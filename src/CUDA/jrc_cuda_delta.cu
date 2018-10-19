/**
 * jrc_cuda_delta.cu
 * block loading delta calculation. should be much faster
 * system('nvcc -ptx -m 64 -arch sm_35 jrc_cuda_rho.cu')
 * iA is multiple of CHUNK (16)
 * J. James Jun, Vidrio Technologies, LLC., 2017 Jun 11
*/

#include <cuda_runtime.h>
// #include "cublas_v2.h"
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define NC 45 // number of Channels
#define CHUNK 16 //previously defined as CHUNK
#define SINGLE_INF (3.402E+38)

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * Step through one B at a time
 * 7/13/17: fDc_spk option added, which uses spike-specific distance cut-off (dc)
 */

// % Matlab syntax
// mrDist12_ = eucl2_dist_(mrFet12, mrFet12(:,1:n1));  %not sqrt
// mlRemove12_ = bsxfun(@ge, viiRho12_ord, viiRho12_ord(1:n1)') ...
//     | abs(bsxfun(@minus, viiSpk12_ord_, viiSpk12_ord_(1:n1)')) > dn_max;
// mrDist12_(mlRemove12_) = nan;
// [vrDelta1, viNneigh1] = min(mrDist12_);

__global__ void jrc_cuda_delta(float * vrDelta1, unsigned int * viNneigh1, const float * mrFet12, const int * viiSpk12_ord, const int * viiRho12_ord, const int * vnConst, const float dc2){
    // int iA = blockIdx.x * CHUNK;    
    int i1 = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK;   // base index of i1
    int tx = threadIdx.x;
    int i1_tx = i1 + tx;
    int n1 = vnConst[0];
    int n12 = vnConst[1];
    int nC = vnConst[2];
    int dn_max = vnConst[3];    
    int fDc_spk = vnConst[4];
    
    __shared__ int viiSpk1_ord_[CHUNK];
    __shared__ int viiRho1_ord_[CHUNK];
    __shared__ float mrFet1_[NC][CHUNK];
    __shared__ float mrDelta1_[NTHREADS][CHUNK];
    __shared__ unsigned int miNneigh1_[NTHREADS][CHUNK]; 
    __shared__ float vrDc1_[CHUNK];  // use if fDc_spk=1
    
    // cache shared memory
    if (tx < nC){ //use tx as iC
        for (int i_c = 0; i_c < CHUNK; ++i_c){
            int i1_c = i_c + i1;
            if (i1_c < n1){
                mrFet1_[tx][i_c] = mrFet12[tx + i1_c * nC];
            }else{
                mrFet1_[tx][i_c] = 0.0f;
            }
        }
    }  
    if (tx < CHUNK && i1_tx < n1){
        viiSpk1_ord_[tx] = viiSpk12_ord[i1_tx];
        viiRho1_ord_[tx] = viiRho12_ord[i1_tx];
    }

    float vr_minDist1[CHUNK];
    unsigned int vi_minIdx1[CHUNK];
    for (int i_c = 0; i_c < CHUNK; ++i_c){
        vr_minDist1[i_c] = SINGLE_INF;
        vi_minIdx1[i_c] = i1 + i_c; // self
    }    
    
    // calculate spike-specific distance cut-off vrDc1_ only if fDc_spk==1
    if (tx < CHUNK && fDc_spk==1){
        vrDc1_[tx] = 0.0f; //init
        //for (int iC = 0; iC < 1; ++iC){ //center only scale
        for (int iC = 0; iC < nC; ++iC){
            float temp_ = mrFet1_[iC][tx];
            vrDc1_[tx] += (temp_ * temp_);
        }
        vrDc1_[tx] *= dc2;
    }
    
    __syncthreads();  
    
    
    // fill in the shared memory A
    for (int i12_tx = tx; i12_tx < n12; i12_tx += blockDim.x){
    //for (int i12_tx = 1; i12_tx < n12; ++i12_tx){
        // compute time difference
        char vlDist_c[CHUNK];
        int iiSpk12_ord_tx = viiSpk12_ord[i12_tx];
        int iiRho12_ord_tx = viiRho12_ord[i12_tx];
        for (int i_c = 0; i_c < CHUNK; ++i_c){
            char di_rho_ = (iiRho12_ord_tx < viiRho1_ord_[i_c]);
            int di_spk_ = ABS(viiSpk1_ord_[i_c] - iiSpk12_ord_tx);
            vlDist_c[i_c] = (di_spk_ <= dn_max) && di_rho_;
        }
        
        // compute distance
        float vrDist_c[CHUNK];
        for (int i_c = 0; i_c < CHUNK; ++i_c) vrDist_c[i_c] = 0.0f;        
        for (int iC = 0; iC < nC; ++iC){
            float fet12_tx = mrFet12[iC + i12_tx * nC];
            for (int i_c = 0; i_c < CHUNK; ++i_c){
                float temp = fet12_tx - mrFet1_[iC][i_c];
                vrDist_c[i_c] += temp * temp;
            }            
        }
        
        // Compare the index and distance
        for (int i_c = 0; i_c < CHUNK; ++i_c){            
            if (vrDist_c[i_c] < vr_minDist1[i_c]){
                if (vlDist_c[i_c] == 1){                
                    vr_minDist1[i_c] = vrDist_c[i_c];
                    vi_minIdx1[i_c] = i12_tx;
                }
            }
        }
    } // while
    
    // collect result from each thread
    for (int i_c = 0; i_c < CHUNK; ++i_c){        
        mrDelta1_[tx][i_c] = vr_minDist1[i_c];
        miNneigh1_[tx][i_c] = vi_minIdx1[i_c];
    }
    __syncthreads();    
    
    // final count    
    //if (tx < CHUNK && i1_tx < n1){
    if (tx < CHUNK){
        float minDist1 = SINGLE_INF;
        unsigned int minIdx1 = i1_tx;
        for (int tx1=0; tx1<blockDim.x; ++tx1){
            if (mrDelta1_[tx1][tx] < minDist1){
                minDist1 = mrDelta1_[tx1][tx];
                minIdx1 = miNneigh1_[tx1][tx];
            }
        }
        //vrDelta1[i1_tx] = sqrtf(minDist1);
        if (i1_tx < n1){
            // vrDelta_ = sqrt(abs(single(vrDelta_) / vrDc2_site(iSite))); %normalize and convert dist
            if (fDc_spk==0){
                vrDelta1[i1_tx] = sqrtf(ABS(minDist1) / dc2); 
            }else{
                vrDelta1[i1_tx] = sqrtf(ABS(minDist1) / vrDc1_[tx]); 
                //vrDelta1[i1_tx] = sqrtf(ABS(minDist1));
            }
            viNneigh1[i1_tx] = minIdx1 + 1; //Matlab index output
        }
    }

} // func