/* ========================================================================== */
/* === GPU/dgemm_custom_simple_1block_batch.cu ============================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/GPU/dgemm_custom_simple_1block_batch.cu.
 * Copyright (C) 2016, Timothy A. Davis
 * CHOLMOD/GPU/dgemm_custom_simple_1block_batch.cu and the CHOLMOD GPU Module 
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


#define  blk_m 16
#define  blk_n 16 
#define  blk_k 16


/* function declarations */
__global__ void dgemm_custom_simple_1block_batch_kernel( int transa, int transb, 
		  				         int *mlist, int *nlist, int *klist,
						         double alpha, 
						         const double **Alist, int *ldalist,				 
					   	         const double **Blist, int *ldblist,
						         double beta,
			  		                 double **Clist, int *ldclist);
			             


extern "C" {
  /* DGEMM wrapper */
void dgemm_custom_simple_1block_batch (cudaStream_t stream, 
  		     	               cublasOperation_t transa, 
				       cublasOperation_t transb, 
				       int *mlist, int *nlist, int *klist, 
				       const double *alpha, 
				       const double **Alist, int *ldalist, 
				       const double **Blist, int *ldblist, 
				       const double *beta, 
				       double **Clist, int *ldclist, int nbatch)                                       
  {

    /* set transpose */
    int ta = (transa == CUBLAS_OP_T) || (transa == CUBLAS_OP_C);
    int tb = (transb == CUBLAS_OP_T) || (transb == CUBLAS_OP_C);

    /* grid size (batch size) */
    dim3 grid(nbatch);   

    /* block size */
    dim3 threads( 16, 16);

    /* batched kernel */
    dgemm_custom_simple_1block_batch_kernel<<< grid, threads, 0, stream >>> (ta, tb, mlist, nlist, klist, *alpha, Alist, ldalist, Blist, ldblist, *beta, Clist, ldclist);
   
  }                                      
}


/* batched DGEMM kernel */
__global__ void dgemm_custom_simple_1block_batch_kernel( int transa, int transb,
		                                         int *mlist, int *nlist, int *klist,
                    			                 double alpha,
		                                         const double **Alist, int *ldalist,
                    			                 const double **Blist, int *ldblist,
		                                         double beta,
                    			                 double **Clist, int *ldclist)
{

    /* allocate shared memory */
    __shared__ double As[blk_m][blk_k+1];
    __shared__ double Bs[blk_k][blk_n+1];

    /* set thread index */
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;

    /* local variables */ 
    int im, in, ik, i, id;
    int colA, colB;
    int rowA, rowB;
    double ABSUM = 0.0;
 
    /* get dimensions for current block (DGEMM) */
    int m = mlist[blockIdx.x];
    int n = nlist[blockIdx.x];
    int k = klist[blockIdx.x];
    int lda = ldalist[blockIdx.x];
    int ldb = ldblist[blockIdx.x];
    int ldc = ldclist[blockIdx.x];
    const double *A = Alist[blockIdx.x];
    const double *B = Blist[blockIdx.x];
    double *C = Clist[blockIdx.x];

    __syncthreads();   

    /* early exit */
    if(!m || !n || !k) return;
   
 
    if (transa) colA = ty;
    else        rowA = tx;

    /* loop over m tiles */
    for(im = 0; im < m; im += blk_m) {

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
            if((rowA < k) && (colA < m)) {        
              As[ty][tx] = A[rowA + lda*colA];  
            }
            else {
              As[ty][tx] = 0.0;
            }
          }
          else {           
            if((rowA < m) && (colA < k)){     
              As[tx][ty] = A[rowA + lda*colA];                 
            }
            else {
              As[tx][ty] = 0.0;
            }   
          }

          /* load B */
          if(transb) {          
            if((rowB < n) && (colB < k)) {
              Bs[ty][tx] = B[rowB + ldb*colB];       
            }
            else {
              Bs[ty][tx] = 0.0;
            }   
          }
          else {                      
            if ((rowB < k) && (colB < n)) {
              Bs[tx][ty] = B[rowB + ldb*colB];          
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
        if (((tx + im) < m) && ((ty + in) < n)) {
          id = tx + im  + ldc*(ty + in);
          C[id] = alpha*ABSUM + beta*C[id];
        }

        if (transb) rowB += blk_n;
        else        colB += blk_n;

      } /* end n tile loop */

      if (transa) colA += blk_m;
      else        rowA += blk_m;

    } /* end m tile loop */

}


