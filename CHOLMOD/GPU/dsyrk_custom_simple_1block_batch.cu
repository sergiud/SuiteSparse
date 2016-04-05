/* ========================================================================== */
/* === GPU/dsyrk_custom_simple_1block_batch.cu ============================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/GPU/dsyrk_custom_simple_1block_batch.cu.
 * Copyright (C) 2016, Timothy A. Davis
 * CHOLMOD/GPU/dsyrk_custom_simple_1block_batch.cu and the CHOLMOD GPU Module 
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

#define  blk_n 16 
#define  blk_k 16


/* function declarations */
__global__ void dsyrk_custom_simple_1block_batch_kernel( int transa,
                                                         int *nlist, int *klist,
                                                         double alpha,
                                                         const double **Alist, int *ldalist,
                                                         double beta,
                                                         double **Clist, int *ldclist);
			             


extern "C" {
  /* DSYRK wrapper */
void dsyrk_custom_simple_1block_batch(cudaStream_t stream,
				      cublasFillMode_t uplo,
                                      cublasOperation_t transa,
                                      int *nlist, int *klist,
                                      const double *alpha,
                                      const double **Alist, int *ldalist,
                                      const double *beta,
                                      double **Clist, int *ldclist, int nbatch)
  {
    /* set transpsoe */
    int  ta = (transa == CUBLAS_OP_T);

    /* check if supported */
    if(uplo != CUBLAS_FILL_MODE_LOWER || transa == CUBLAS_OP_C ) {
      printf("only the following properties are supported:\n");
      printf("uplo  = CUBLAS_FILL_MODE_LOWER\n");
      printf("trans = CUBLAS_OP_N, CUBLAS_OP_T\n");
    }

    /* grid size (batch size) */
    dim3 grid( nbatch );

    /* block size */
    dim3 threads( blk_n, blk_n);
 
    /* batched kernel */
    dsyrk_custom_simple_1block_batch_kernel<<< grid, threads, 0, stream >>> (ta, nlist, klist, *alpha, Alist, ldalist, *beta, Clist, ldclist);

  }
                                 
}


/* batched DSYRK kernel */
__global__ void dsyrk_custom_simple_1block_batch_kernel( int transa,
                                                         int *nlist, int *klist,
                                                         double alpha,
                                                         const double **Alist, int *ldalist,
                                                         double beta,
                                                         double **Clist, int *ldclist)
{

    /* allocate shared memory */
    __shared__ double As[blk_n][blk_k+1];
    __shared__ double Bs[blk_k][blk_n+1];

    /* set thread index */
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;

    /* local variables */ 
    int im, in, ik, i, id;
    int colA, colB;
    int rowA, rowB;
    double ABSUM = 0.0;

    int transb = !transa;
 
    /* get dimensions for current block (DSYRK) */
    int n = nlist[blockIdx.x];
    int k = klist[blockIdx.x];
    int lda = ldalist[blockIdx.x];
    int ldc = ldclist[blockIdx.x];
    const double *A = Alist[blockIdx.x];
    double *C = Clist[blockIdx.x];

    __syncthreads();

    /* early exit */
    if(!n || !k) return;

 
    if (transa) colA = ty;
    else        rowA = tx;

    /* loop over m tiles */
    for(im = 0; im < n; im += blk_n) {

      if (transb) rowB =  tx;
      else colB =  ty;

      /* loop over n tiles */
      for(in = 0; in < n; in += blk_n) {

        ABSUM = 0.0;

        if (transa) rowA = tx;
        else 	    colA = ty;

        if (transb) colB =  ty;
        else 	    rowB =  tx;

        /* synchronize threads */
        __syncthreads();
    
        /* loop over k tiles */
        for(ik = 0; ik < k; ik += blk_k) { 

          /* load A */
          if(transa) {
            if((rowA < k) && (colA < n)) {        
              As[ty][tx] = A[rowA + lda*colA];  
            }
            else {
              As[ty][tx] = 0.0;
            }
          }
          else {           
            if((rowA < n) && (colA < k)) {     
              As[tx][ty] = A[rowA + lda*colA];                 
            }
            else {
              As[tx][ty] = 0.0;
            }   
          }

          /* load B */
          if(transb) {          
            if((rowB < n) && (colB < k)) {
              Bs[ty][tx] = A[rowB + lda*colB];       
            }
            else {
              Bs[ty][tx] = 0.0;
            }   
          }
          else {                      
            if ((rowB < k) && (colB < n)) {
              Bs[tx][ty] = A[rowB + lda*colB];          
            }
            else {
              Bs[tx][ty] = 0.0;
            }
          } 

          /* synchronize threads */
          __syncthreads();

          /* compute A*B */
          for(i = 0; i < blk_k; i++){
            ABSUM += As[tx][i]*Bs[i][ty];
          }

          if (transa) rowA += blk_k;
          else 	      colA += blk_k;
          if (transb) colB += blk_k;
          else 	      rowB += blk_k;

          /* synchronize threads */
          __syncthreads();

        }  /* end k tile loop */
        

        /* compute C = alpha*A*B + beta*C */
        if (((tx + im) < n) && ((ty + in) <= (tx + im))) {
          id = tx + im  + ldc*(ty + in);
          C[id] = alpha*ABSUM + beta*C[id];
        }

        if (transb) rowB += blk_n;
        else        colB += blk_n;

      } /* end n tile loop */

      if (transa) colA += blk_n;
      else        rowA += blk_n;

    } /* end m tile loop */

}


