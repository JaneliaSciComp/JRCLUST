/**
 * jrc_cuda_rho.cu
 * block loading rho calculation. should be much faster
 * system('nvcc -ptx -m 64 -arch sm_35 jrc_cuda_rho.cu')
 * i1 is multiple of chunk (16)
 * J. James Jun, Vidrio Technologies, LLC., 2017 Jun 11
 * 7/13/17: fDc_spk option added, which uses spike-specific distance cut-off (dc)
*/

#include <cuda_runtime.h>
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define MAXDIM 45
#define CHUNKSIZE 16
#define SINGLE_INF (3.402E+38) // equivalent to NAN. consider -1 value

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * Step through one B at a time
 */
__global__ void jrc_cuda_rho(float *rho, const float *site_features, const int *spike_order, const int *site_constants, const float dist_cut2) {
    int i1 = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNKSIZE; // base index of i1
    int thread_x = threadIdx.x; // nThreadsGPU for i12 index
    int i1_thread_x = i1 + thread_x;
    int n_spikes_primary = site_constants[0];
    int n_spikes_all = site_constants[1];
    int n_features = site_constants[2];
    int dn_max = site_constants[3];
    int fDc_spk = site_constants[4];

    __shared__ int spike_order_chunk[CHUNKSIZE];
    __shared__ float features_primary[MAXDIM][CHUNKSIZE];
    __shared__ int rho_chunk[NTHREADS][CHUNKSIZE]; // count then divide later
    __shared__ int mnComp1_[NTHREADS][CHUNKSIZE]; // count number of elements compared
    __shared__ float vrDc1_[CHUNKSIZE];  // use if fDc_spk=1

    // cache shared memory
    if (thread_x < n_features) {
        for (int i_c = 0; i_c < CHUNKSIZE; i_c++) {
            int i1_c = i_c + i1;
            if (i1_c < n_spikes_primary) {
                features_primary[thread_x][i_c] = site_features[thread_x + i1_c * n_features];
            } else {
                features_primary[thread_x][i_c] = 0.0f;
            }
        }
    }

    if (thread_x < CHUNKSIZE && i1_thread_x < n_spikes_primary) {
        spike_order_chunk[thread_x] = spike_order[i1_thread_x];
    }

    // initialize rho
    for (int i_c = 0; i_c < CHUNKSIZE; i_c++) {
        rho_chunk[thread_x][i_c] = 0;
        mnComp1_[thread_x][i_c] = 0;
    }

    // calculate spike-specific distance cut-off vrDc1_ only if fDc_spk==1
    if (thread_x < CHUNKSIZE && fDc_spk == 1) {
        vrDc1_[thread_x] = 0.0f; //init
        for (int i_feature = 0; i_feature < n_features; i_feature++) {
            float temp = features_primary[i_feature][thread_x];
            vrDc1_[thread_x] += (temp * temp);
        }
        vrDc1_[thread_x] *= dist_cut2;
    }

    __syncthreads();

    // Inspect distance relationship between i1 and i12_tx
    for (int i12_tx = thread_x; i12_tx < n_spikes_all; i12_tx += blockDim.x) {
        int iiSpk12_ord_tx = spike_order[i12_tx];

        // compute distance
        float feature_dists2_chunk[CHUNKSIZE]; // square of pairwise feature distances for chunk
        for (int i_c = 0; i_c < CHUNKSIZE; i_c++) {
            feature_dists2_chunk[i_c] = 0.0f;
        }

        for (int i_feature = 0; i_feature < n_features; i_feature++) {
            float fet12_tx = site_features[i_feature + i12_tx * n_features];
            for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
                float temp = fet12_tx - features_primary[i_feature][i_c]; // z_i = x_i - y_i
                feature_dists2_chunk[i_c] += temp * temp;                 // dist += z_i^2
            }
        }

        // Compare the index and distance
        for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
            int time_dist = ABS(spike_order_chunk[i_c] - iiSpk12_ord_tx);
            if (time_dist <= dn_max) {
                ++mnComp1_[thread_x][i_c];
                if (fDc_spk == 0) {
                    if (feature_dists2_chunk[i_c] <= dist_cut2) {
                        ++rho_chunk[thread_x][i_c];
                    }
                } else {
                    if (feature_dists2_chunk[i_c] < vrDc1_[i_c]) {
                        ++rho_chunk[thread_x][i_c];
                    }
                }
            }
        }
    } // for

    // final count
    __syncthreads();

    if (thread_x < CHUNKSIZE) {  // use thread_x as i_c
        int nRho1 = 0;
        int nComp1 = 0;
        for (int tx1 = 0; tx1 < blockDim.x; tx1++) {
            nRho1 += rho_chunk[tx1][thread_x];
            nComp1 += mnComp1_[tx1][thread_x];
        }

        if (i1_thread_x < n_spikes_primary) {
            rho[i1_thread_x] = (float)(((double) (nRho1)) / ((double) nComp1));
        }
    }
}