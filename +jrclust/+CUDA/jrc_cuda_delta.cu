/**
 * jrc_cuda_delta.cu
 * block loading delta calculation. should be much faster
 * system('nvcc -ptx -m 64 -arch sm_35 jrc_cuda_rho.cu')
 * iA is multiple of CHUNKSIZE (16)
 * J. James Jun, Vidrio Technologies, LLC., 2017 Jun 11
*/

#include <cuda_runtime.h>
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define MAXDIM 45 // number of Channels
#define CHUNKSIZE 16 // previously defined as CHUNK
#define SINGLE_INF (3.402E+38)

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * Step through one B at a time
 * 7/13/17: fDc_spk option added, which uses spike-specific distance cut-off (dc)
 */
__global__ void jrc_cuda_delta(float *delta, unsigned int *nneigh, const float *site_features, const int *spike_order, const int *rho_order, const int *site_constants, const float dist_cut2) {
    int i1 = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNKSIZE; // base index of i1
    int thread_x = threadIdx.x;
    int i1_thread_x = i1 + thread_x;
    int n_spikes_primary = site_constants[0];
    int n_spikes_all = site_constants[1];
    int n_features = site_constants[2];
    int time_dist_cut = site_constants[3];
    int fDc_spk = site_constants[4];

    __shared__ int spike_order_chunk[CHUNKSIZE];
    __shared__ int rho_order_chunk[CHUNKSIZE];
    __shared__ float features_primary[MAXDIM][CHUNKSIZE];
    __shared__ float mrDelta1_[NTHREADS][CHUNKSIZE];
    __shared__ unsigned int miNneigh1_[NTHREADS][CHUNKSIZE];
    __shared__ float vrDc1_[CHUNKSIZE]; // use if fDc_spk == 1

    // cache shared memory, 1/2
    if (thread_x < n_features) { // use thread_x as iC
        for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
            int i1_c = i_c + i1;
            if (i1_c < n_spikes_primary) {
                features_primary[thread_x][i_c] = site_features[thread_x + i1_c * n_features];
            } else {
                features_primary[thread_x][i_c] = 0.0f;
            }
        }
    }

    // cache shared memory, 2/2
    if (thread_x < CHUNKSIZE && i1_thread_x < n_spikes_primary) {
        spike_order_chunk[thread_x] = spike_order[i1_thread_x];
        rho_order_chunk[thread_x] = rho_order[i1_thread_x];
    }

    float mindist_chunk[CHUNKSIZE];
    unsigned int nneigh_chunk[CHUNKSIZE];
    for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
        mindist_chunk[i_c] = SINGLE_INF;
        nneigh_chunk[i_c] = i1 + i_c; // self
    }

    // calculate spike-specific distance cut-off vrDc1_ only if fDc_spk == 1
    if (thread_x < CHUNKSIZE && fDc_spk == 1) {
        vrDc1_[thread_x] = 0.0f; // init
        for (int iC = 0; iC < n_features; ++iC) {
            float temp_ = features_primary[iC][thread_x];
            vrDc1_[thread_x] += (temp_ * temp_);
        }
        vrDc1_[thread_x] *= dist_cut2;
    }

    __syncthreads();

    // fill in the shared memory A
    for (int i12_tx = thread_x; i12_tx < n_spikes_all; i12_tx += blockDim.x) {
        // compute time difference
        char nearby_in_time[CHUNKSIZE];
        int i_spike_order = spike_order[i12_tx];
        int i_rho_order = rho_order[i12_tx];

        for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
            char rho_is_larger = (i_rho_order < rho_order_chunk[i_c]); // is rho larger?
            int time_dist = ABS(spike_order_chunk[i_c] - i_spike_order); // is the spike nearby in time?
            nearby_in_time[i_c] = (time_dist <= time_dist_cut) && rho_is_larger;
        }

        // compute distance
        float feature_dists2_chunk[CHUNKSIZE]; // square of pairwise feature distances for chunk
        for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
            feature_dists2_chunk[i_c] = 0.0f;
        }

        for (int iC = 0; iC < n_features; ++iC) {
            float fet12_tx = site_features[iC + i12_tx * n_features];
            for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
                float temp = fet12_tx - features_primary[iC][i_c]; // z_i = x_i - y_i
                feature_dists2_chunk[i_c] += temp * temp;          // dist += z_i^2
            }
        }

        // Compare the index and distance
        for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
            if ((nearby_in_time[i_c] == 1) && (feature_dists2_chunk[i_c] < mindist_chunk[i_c])) {
                mindist_chunk[i_c] = feature_dists2_chunk[i_c];
                nneigh_chunk[i_c] = i12_tx;
            }
        }
    } // for

    // collect result from each thread
    for (int i_c = 0; i_c < CHUNKSIZE; ++i_c) {
        mrDelta1_[thread_x][i_c] = mindist_chunk[i_c];
        miNneigh1_[thread_x][i_c] = nneigh_chunk[i_c];
    }
    __syncthreads();

    // final count
    if (thread_x < CHUNKSIZE) {
        float minDist1 = SINGLE_INF;
        unsigned int minIdx1 = i1_thread_x;
        for (int tx1 = 0; tx1 < blockDim.x; ++tx1) {
            if (mrDelta1_[tx1][thread_x] < minDist1) {
                minDist1 = mrDelta1_[tx1][thread_x];
                minIdx1 = miNneigh1_[tx1][thread_x];
            }
        }

        if (i1_thread_x < n_spikes_primary) {
            if (fDc_spk == 0) {
                delta[i1_thread_x] = sqrtf(ABS(minDist1) / dist_cut2);
            } else {
                delta[i1_thread_x] = sqrtf(ABS(minDist1) / vrDc1_[thread_x]);
            }
            nneigh[i1_thread_x] = minIdx1 + 1; // Matlab index output
        }
    }
} // func