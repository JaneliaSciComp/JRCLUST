/**
 * block loading rho calculation. should be much faster
 * system('nvcc -ptx citydist_rho4.cu')
 * iA is multiple of CHUNK (16)
*/

#include <cuda_runtime.h>
// #include "cublas_v2.h"
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define NC (1+6*2) // number of Channels
#define CHUNK 16 //previously defined as CHUNK
#define SINGLE_INF (3.402E+38)

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * Step through one B at a time
 */
__global__ void citydist_sorted_delta(float const *A, unsigned int const *I, float *D, unsigned int *N, int const nA, int const nneigh, int const nC){
    // int iA = blockIdx.x * CHUNK;    
    int iA = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK;
    int iA1;
    int tx = threadIdx.x;
    float vr_minDist1[CHUNK];
    unsigned int vi_minIdx1[CHUNK];
    __shared__ unsigned int svi_I_A1[CHUNK];
    __shared__ float smr_A1[NC][CHUNK];
    __shared__ float smr_delta1[NTHREADS][CHUNK];
    __shared__ unsigned int smi_nneigh1[NTHREADS][CHUNK]; 
    
    // cache A
    int iC = tx;    
    if (iC < nC){ //use tx as iC
        for (iA1 = 0; iA1 < CHUNK; ++iA1){
            if (iA + iA1 < nA){
                smr_A1[iC][iA1] = A[iC + (iA+iA1)*nC];
            }else{
                smr_A1[iC][iA1] = SINGLE_INF;
            }   
        }
    }
    
    iA1 = tx; // batch index
    if (iA1 < CHUNK){
        if (iA + iA1 < nA){
            svi_I_A1[iA1] = I[iA + iA1];  
        }else{
            svi_I_A1[iA1] = nA + 1;  // out of range
        }
    }
    for (iA1 = 0; iA1 < CHUNK; ++iA1){
        vr_minDist1[iA1] = SINGLE_INF;
        vi_minIdx1[iA1] = iA + iA1;
    }    
    __syncthreads();

    // fill in the shared memory A
    int iB_min = MAX(iA - nneigh, 0);
    int iB_max = MIN(iA + nneigh + CHUNK - 1, nA-1);
    if (nneigh==0){
        iB_min = 0; 
        iB_max = nA-1;
    }      
    int iB = iB_min + tx;  
    while (iB <= iB_max){
        float vr_dist1[CHUNK];
        for (iA1 = 0; iA1 < CHUNK; ++iA1) vr_dist1[iA1] = 0.0f;
        for (iC = 0; iC < nC; ++iC){
            float Btemp = A[iC + iB*nC];
            for (iA1 = 0; iA1 < CHUNK; ++iA1){
                float temp = Btemp - smr_A1[iC][iA1];
                vr_dist1[iA1] += ABS(temp);
            }            
        }          
        unsigned int IiB = I[iB];
        for (iA1 = 0; iA1 < CHUNK; ++iA1){            
            if (vr_dist1[iA1] < vr_minDist1[iA1]){
                if (IiB < svi_I_A1[iA1]){                
            //if (vr_dist1[iA1] < vr_minDist1[iA1] && vr_dist1[iA1]>0){
                //if (IiB < svi_I_A1[iA1] && iB != iA+iA1){
                    int dab = ABS(iA + iA1 - iB);
                    if (dab <= nneigh || nneigh==0){
                        vr_minDist1[iA1] = vr_dist1[iA1];
                        vi_minIdx1[iA1] = iB;
                    }
                }
            }
        }
        iB += blockDim.x;
    } // while
    
    // collect result from each thread
    for (iA1 = 0; iA1 < CHUNK; ++iA1){        
        smr_delta1[tx][iA1] = vr_minDist1[iA1];
        smi_nneigh1[tx][iA1] = vi_minIdx1[iA1];
    }
    __syncthreads();    
    
    // final count    
    iA1 = tx;
    if (iA1 < CHUNK && iA + iA1 < nA){
        float minDist1 = SINGLE_INF;
        unsigned int minIdx1 = iA + iA1;
        for (int tx1=0; tx1<blockDim.x; ++tx1){
            if (smr_delta1[tx1][iA1] < minDist1){
                minDist1 = smr_delta1[tx1][iA1];
                minIdx1 = smi_nneigh1[tx1][iA1];
            }
        }
        D[iA + iA1] = sqrtf(minDist1);
        N[iA + iA1] = minIdx1;        
    }

} // func