/* ========================================================================== */
/* === GPU/dtrsm_custom_simple_1block_batch.cu ============================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/GPU/dtrsm_custom_simple_1block_batch.cu.
 * Copyright (C) 2016, Timothy A. Davis
 * CHOLMOD/GPU/dtrsm_custom_simple_1block_batch.cu and the CHOLMOD GPU Module 
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


#define  blk_m 64	/* # threads in block */

#define A(I,J) A[(I) + ((J))* (lda)]
#define B(I,J) B[(I) + ((J))* (ldb)]

/* function declarations */
__global__ void dtrsm_lower_right_1block_batch_kernel( int trans, 
		  		    	               int *mlist, int *nlist,
				                       double alpha, 
				    	               const double **Alist, int *ldalist,				 
 			    	                       double **Blist, int *ldblist);
			  
__global__ void dtrsm_lower_right_transpose_1block_batch_kernel( int trans,
                                                                 int *mlist, int *nlist,
                                                                 double alpha,
                                                                 const double **Alist, int *ldalist,
                                                                 double **Blist, int *ldblist);





extern "C" {
  /* DTRSM wrapper */
void dtrsm_custom_simple_1block_batch (cudaStream_t stream,
				       cublasSideMode_t side, 
                                       cublasFillMode_t uplo, 
                                       cublasOperation_t trans, 
                                       cublasDiagType_t diag, 
                                       int *mlist, int *nlist,
		    	               const double *alpha, 
                                       const double **Alist, int *ldalist, 
			               double **Blist, int *ldblist, int nbatch)                                       
  {
    /* set transpose */
    int ta = (trans == CUBLAS_OP_T);

    /* check if supported */
    if(side != CUBLAS_SIDE_RIGHT || uplo != CUBLAS_FILL_MODE_LOWER || diag != CUBLAS_DIAG_NON_UNIT) {
      printf("only the following properties are supported:\n");
      printf("side  = CUBLAS_SIDE_RIGHT\n");
      printf("uplo  = CUBLAS_FILL_MODE_LOWER\n");
      printf("trans = CUBLAS_OP_N, CUBLAS_OP_T\n");
      printf("diag  = CUBLAS_DIAG_NON_UNIT\n");
    }  

    /* grid size (batch size) */
    dim3 grid(nbatch);

    /* block size */
    dim3 threads(blk_m);

    /* batched kernel */
    if(!ta) dtrsm_lower_right_1block_batch_kernel<<< grid, threads, 0, stream>>> (ta, mlist, nlist, *alpha, Alist, ldalist, Blist, ldblist);
    else    dtrsm_lower_right_transpose_1block_batch_kernel<<< grid, threads, 0, stream>>> (ta, mlist, nlist, *alpha, Alist, ldalist, Blist, ldblist);
   
  }                                      
}

/* batched DTRSM kernel (non transpose)  */
__global__ void  dtrsm_lower_right_1block_batch_kernel( int trans,
		                                        int *mlist, int *nlist,
                    		                        double alpha,
		                                        const double **Alist, int *ldalist,
                    		                        double **Blist, int *ldblist)
{
  /* set thread index */
  const int tx = threadIdx.x;

  /* local variables */
  int j,k,im;
  double temp;
  
  /* get dimensions of current block (DTRSM) */
  int m = mlist[blockIdx.x];
  int n = nlist[blockIdx.x];
  int lda = ldalist[blockIdx.x];
  int ldb = ldblist[blockIdx.x];
  const double *A = Alist[blockIdx.x];
  double *B = Blist[blockIdx.x];

  __syncthreads();

  /* early exit */
  if(!m || !n) return; 


  /* loop over m tiles */
  for(im = 0; im < m; im += blk_m) {

    if((tx + im) < m) {
      for (j = n - 1; j >= 0; j--) {

        /* compute B = B - sum(A*B) */
        for (k = j+1; k < n; k++) {
          temp = A(k,j);
          B(tx+im,j) -= temp * B(tx+im,k);
        }

        /* non unit */      
        temp = 1.0/A(j,j);

        /* compute B = alpha*B/A */
        B(tx+im,j) = alpha * temp * B(tx+im,j);
      }
    }
  
    /* syncronize threads */
    __syncthreads();

  } /* end m tile loop */


}






/* batched DTRSM kernel (transpose) */
__global__ void  dtrsm_lower_right_transpose_1block_batch_kernel( int trans,
                                                                  int *mlist, int *nlist,
                                                                  double alpha,
                                                                  const double **Alist, int *ldalist,
                                                                  double **Blist, int *ldblist)
{
  /* thread index */
  const int tx = threadIdx.x;

  /* local variables */
  int j,k,im;
  double temp;

  /* get dimensions of current block (DTRSM) */
  int m = mlist[blockIdx.x];
  int n = nlist[blockIdx.x];
  int lda = ldalist[blockIdx.x];
  int ldb = ldblist[blockIdx.x];
  const double *A = Alist[blockIdx.x];
  double *B = Blist[blockIdx.x];

  __syncthreads();


  /* early exit */
  if(!m || !n) return;


  /* loop over m tiles */
  for(im = 0; im < m; im += blk_m) {

    if((tx + im) < m) {
      for (k = 0; k < n; k++) {

        /* non unit */
        temp = 1.0/A(k,k);

        /* compute B = alpha*B/A */
        B(tx+im,k) = alpha * temp * B(tx+im,k);

        /* compute B = B - sum(A*B) */
        for (j = k + 1; j < n; j++) {
          temp = A(j,k);
          B(tx+im,j) -= temp * B(tx+im,k);
        }

      }
    } 

    /* syncronize threads */
    __syncthreads();

  } /* end m tiles loop */

}
