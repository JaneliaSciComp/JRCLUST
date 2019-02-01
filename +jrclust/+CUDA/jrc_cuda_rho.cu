/**
 * jrc_cuda_rho.cu
 * block loading rho calculation. should be much faster
 * system('nvcc -ptx -m 64 -arch sm_35 jrc_cuda_rho.cu')
 * i1 is multiple of chunk (16)
 * J. James Jun, Vidrio Technologies, LLC., 2017 Jun 11
 * 7/13/17: fDc_spk option added, which uses spike-specific distance cut-off (dc)
*/

#include <cuda_runtime.h>
// #include "cublas_v2.h"
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define NC 45 //max dimm
#define CHUNK 16
#define SINGLE_INF (3.402E+38) // equipvalent to NAN. consider -1 value

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * Step through one B at a time
 */
__global__ void jrc_cuda_rho(float * vrRho1, const float * mrFet12, const int * viiSpk12_ord, const int *  vnConst, const float dc2){
//__global__ void jrc_cuda_rho(int *vnRho1, int *vnComp1, float const *mrFet12, int const *viiSpk12_ord, int const *vnC4, float const dc2){
    int i1 = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK;   // base index of i1    
    int tx = threadIdx.x; //nThreadsGPU for i12 index    
    int i1_tx = i1+tx;
    int n1 = vnConst[0];
    int n12 = vnConst[1];
    int nC = vnConst[2];
    int dn_max = vnConst[3];    
    int fDc_spk = vnConst[4];
    
    __shared__ int viiSpk1_ord_[CHUNK];
    __shared__ float mrFet1_[NC][CHUNK];    
    __shared__ int mnRho1_[NTHREADS][CHUNK]; // count then divide later    
    __shared__ int mnComp1_[NTHREADS][CHUNK]; // count number of elements compared        
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
    if (tx < CHUNK && i1_tx < n1) viiSpk1_ord_[tx] = viiSpk12_ord[i1_tx];
    
    for (int i_c = 0; i_c < CHUNK; ++i_c){
        mnRho1_[tx][i_c] = 0; // initialize rho
        mnComp1_[tx][i_c] = 0;
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

    
    // Inspect distance relationship between i1 and i12_tx
    for (int i12_tx = tx; i12_tx < n12; i12_tx += blockDim.x){
    //for (int i12_tx = 1; i12_tx < n12; ++i12_tx){
        // compute time difference
        //char vlDist_c[CHUNK];
        int iiSpk12_ord_tx = viiSpk12_ord[i12_tx];        
        /*for (int i_c = 0; i_c < CHUNK; ++i_c){
            int di_spk_tx = ABS(viiSpk1_ord_[i_c] - iiSpk12_ord_tx);
            vlDist_c[i_c] = (di_spk_tx <= dn_max);
        } */
        
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
            int di_spk_tx = ABS(viiSpk1_ord_[i_c] - iiSpk12_ord_tx);
            if (di_spk_tx <= dn_max){
            //if (vlDist_c[i_c] == 1){
                ++mnComp1_[tx][i_c];
                if (fDc_spk==0){
                    if (vrDist_c[i_c] <= dc2) ++mnRho1_[tx][i_c];
                }else{
                    if (vrDist_c[i_c] < vrDc1_[i_c]) ++mnRho1_[tx][i_c];
                }
            }
        }
    } // while
    
    // final count
    __syncthreads();
    //if (tx < CHUNK && i1_tx < n1){  // use tx as i_c
    if (tx < CHUNK){  // use tx as i_c
        int nRho1 = 0;
        int nComp1 = 0;
        for (int tx1=0; tx1<blockDim.x; ++tx1){
            nRho1 += mnRho1_[tx1][tx];
            nComp1 += mnComp1_[tx1][tx];
        }
        if (i1_tx < n1){
            //if (nRho1<1) nRho1 = 1; 
            vrRho1[i1_tx] = (float)(((double)(nRho1)) / ((double)nComp1));
        }
        // vnRho1[i1 + i_c_] = nRho1 - 1;
        // vnComp1[i1 + i_c_] = nComp1;
    }
    //vnRho1[0] = blockDim.x; //debug
    //vnComp1[0] = blockDim.x; //debug
} // func