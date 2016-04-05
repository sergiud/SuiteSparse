/* ========================================================================== */
/* === GPU/dpotrf_custom_simple_1block_batch.cu ============================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/GPU/dpotrf_custom_simple_1block_batch.cu.
 * Copyright (C) 2016, Timothy A. Davis
 * CHOLMOD/GPU/dpotrf_custom_simple_1block_batch.cu and the CHOLMOD GPU Module 
 * are licensed under Version 2.0 of the GNU General Public License.  See 
 * gpl.txt for a text of the license.  CHOLMOD is also available under other 
 * licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>


#define  blk_n 32	/* # threads in block */

#define A(I,J) A[(I) + ((J))* (lda)]
#define B(I,J) B[(I) + ((J))* (ldb)]

/* function declarations */
__global__ void dpotf2_custom_simple_1block_batch_column_kernel( int *nlist, double **alist, int *ldalist);

			  

extern "C" {


/* DPOTRF wrapper */
void dpotrf_custom_simple_1block_batch(cudaStream_t stream,
				       cublasFillMode_t uplo,
                                       int *nlist, 
                                       double **alist,
	  		               int *ldalist, int *info, int nbatch)                                       
  {

    /* check if supported */
    if(uplo != CUBLAS_FILL_MODE_LOWER) {
      printf("only the following properties are supported:\n");
      printf("uplo  = CUBLAS_FILL_MODE_LOWER\n");
    }  

    /* grid size (batch size) */
    dim3 grid(nbatch);

    /* block size */
    dim3 threads(blk_n);

    /* batched kernel */
    dpotf2_custom_simple_1block_batch_column_kernel<<< grid, threads, 0, stream>>> (nlist, alist, ldalist);
   
  }                                      
}



/* Column Cholesky Method:
 * compute the Cholesky factorization A = L*L' by computing
 * each column using previous columns. (column-major order) */

/* batched dpotrf kernel */
__global__ void dpotf2_custom_simple_1block_batch_column_kernel( int *nlist, double **alist, int *ldalist)
{
  /* thread index */
  const int tx = threadIdx.x;

  /* local variables */
  int i,j,k,in,btx;
  double ajj;

  /* set batch dimensions */
  int n     = nlist[blockIdx.x];
  int lda   = ldalist[blockIdx.x];
  double *a = alist[blockIdx.x];

  __syncthreads();

  /* early exit */
  if(!n) return;

  /* loop over columns */
  for(j = 0; j < n; j++) {

    /* compute l_ij */
    for(in = 0; in < n; in += blk_n) {
      btx = tx + in;
      if(btx > (j - 1) && btx < n) {
        i = btx;
        for(k = 0; k < j; k++) {
          a[i + j*lda] -= a[i + k*lda] * a[j + k*lda];   /* put [a[j + k*lda]] in shared memory.. */
        }
      }
    }
    __syncthreads();

    /* compute l_jj */
    ajj = sqrt( a[j + j*lda] );

    /* compute l_ij */
    for(in = 0; in < n; in += blk_n) {
      btx = tx + in;
      if(btx > j && btx < n) {
        i = btx;
        a[i + j*lda] /= ajj;
      }
    }
    __syncthreads();

    /* store l_jj */
    if(!tx) {
      a[j + j*lda] = ajj;
    }
    __syncthreads();


  }


}




