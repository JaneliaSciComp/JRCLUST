/**
 * block loading rho calculation. should be much faster
 * system('nvcc -ptx citydist_rho4.cu')
 * iA is multiple of chunk (16)
*/

#include <cuda_runtime.h>
// #include "cublas_v2.h"
#include <math.h>
#define ABS(my_val) ((my_val) < 0) ? (-1*(my_val)) : (my_val)
#define MIN(A,B) ((A)<(B)) ? (A) : (B)
#define MAX(A,B) ((A)>(B)) ? (A) : (B)
#define NTHREADS 128
#define NC (1+6*2)
// #define NC (9)
#define CHUNK 16
#define SINGLE_INF (3.402E+38)

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 * Step through one B at a time
 */
__global__ void eucldist_sorted_rho_exp(float const *A, float *D, int const nA, int const nneigh, int const nC, float const dc){
    int iA = (blockIdx.x + blockIdx.y * gridDim.x) * CHUNK;    
    int tx = threadIdx.x;
    __shared__ float sA[NC][CHUNK];
    __shared__ float rho1[NTHREADS][CHUNK];
    
    // cache A
    if (tx < nC){ //use tx as iC
        for (int i=0; i<CHUNK; ++i){
            if (iA+i < nA){
                sA[tx][i] = A[tx + (iA+i)*nC];
            }else{
                sA[tx][i] = SINGLE_INF;
            }   
        }
    }
    for (int i=0; i<CHUNK; ++i) rho1[tx][i] = 0.0f;
    __syncthreads();
    
    // fill in the shared memory A
    float dc2 = dc*dc;
    int iB_min = MAX(iA - nneigh, 0);
    int iB_max = MIN(iA + nneigh + CHUNK - 1, nA-1);
    int iB = iB_min + tx; //MAX(tx, iB_min); // tx is index for B    
    while (iB <= iB_max){
        float dist[CHUNK];
        // calculate distance to B        
        for (int i=0; i<CHUNK; ++i) dist[i] = 0.0f;
        for (int iC=0; iC<nC; ++iC){
            float Btemp = A[iC + iB*nC];
            for (int i=0; i<CHUNK; ++i){
                float temp = Btemp - sA[iC][i];
                dist[i] += temp * temp;
            }            
        }          
        for (int i=0; i<CHUNK; ++i){            
            int dab = ABS(iA+i-iB);
            if (dab<=nneigh){
                if (iA+i < nA && iA+i != iB){
                    rho1[tx][i] += expf(-1*dist[i]/dc2);
                }
            }
        }
        iB += blockDim.x;
    } // while
    
    // final count
    __syncthreads();    
    // if (tx < CHUNK) D[iA+tx] = rho1[tx];
    if (tx < CHUNK){
        float sum = 0.0f;
        for (int tx1=0; tx1<blockDim.x; ++tx1)
            sum += rho1[tx1][tx];
        if (iA+tx<nA) D[iA+tx] = sum;
    }
} // func