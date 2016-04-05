/* ========================================================================== */
/* === GPU/t_cholmod_subtree ================================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/GPU Module.  Copyright (C) 2005-2012, Timothy A. Davis
 * The CHOLMOD/GPU Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/*
 * File:
 *   t_cholmod_subtree
 *
 * Description:
 *   Contains functions for GPU subtrees
 *   algorithm.
 *
 */


#ifdef SUITESPARSE_CUDA


/* include */
#include <string.h>
#include "cholmod_template.h"
#include "batchKernels.h"
#include <cusolverDn.h>
#include <time.h>











/*
 * Function: 
 *   gpu_init
 *
 * Description:
 *   Performs required initialization for GPU computing. 
 *   Returns 0 if there is an error, disabling GPU computing (useGPU = 0)
 *
 */
int TEMPLATE2 (CHOLMOD (gpu_init))
  (
   cholmod_common *Common,
   cholmod_factor *L,
   cholmod_sparse *A,
   cholmod_global_pointers *gb_p,
   cholmod_gpu_pointers *gpu_p,
   int gpuid
   )
{
  /* local variables */
  Int i, k, LsDim, ApDim, AiDim, AxDim;
  Int *h_Ls, *h_Ap, *h_Ai, *Ls, *Ap, *Ai;
  double *h_Ax, *Ax;
  int numThreads;

  cudaError_t cudaErr;
  cublasStatus_t cublasError;
  cusolverStatus_t cusolverErr;


  numThreads = Common->ompNumThreads;

  {

#ifdef USE_NVTX
    nvtxRangeId_t range1 = nvtxRangeStartA("gpu_init");
#endif

    /* set gpu memory pointers */
    gpu_p->gpuPtr[gpuid]		 = Common->dev_mempool[gpuid];

    /* type (double *) */
    
    /* pointers to supernode matrices */
    gpu_p->d_ptrSuper[gpuid] 	 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= 3*gb_p->ptrSuperSize;
    
    /* pointers to descendant matrices */
    gpu_p->d_ptrDesc[gpuid] 	 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= 6*gb_p->ptrDescSize;
    
    /* type (double) */

    /* Lx */
    gpu_p->d_Lx[gpuid] 		 = gpu_p->gpuPtr[gpuid];		
    gpu_p->gpuPtr[gpuid] 		+= gb_p->LxSize;
    
    /* C */
    gpu_p->d_C[gpuid] 		 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= gb_p->CSize;
    
    /* Ax */ 
    gpu_p->d_Ax[gpuid] 		 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= gb_p->AxSize;  
    
    /* type (Int) */
    
    /* Ap */
    gpu_p->d_Ap[gpuid] 		 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= gb_p->ApSize;

    /* Ai */
    gpu_p->d_Ai[gpuid] 		 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= gb_p->AiSize;
    
    /* Ls */
    gpu_p->d_Ls[gpuid] 		 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= gb_p->LsSize;
    
    /* Map */
    gpu_p->d_Map[gpuid] 		 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= gb_p->MapSize;

    /* type (int) */
    
    /* dimensions of supernode matrices */
    gpu_p->d_dimSuper[gpuid] 	 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= 13*gb_p->dimSuperSize;
    
    /* dimensions of descendant matrices */ 
    gpu_p->d_dimDesc[gpuid] 	 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= 14*gb_p->dimDescSize;

    /* info */
    gpu_p->d_info[gpuid] 		 = gpu_p->gpuPtr[gpuid];
    gpu_p->gpuPtr[gpuid] 		+= sizeof(int);

    /* devSync memory for custom cusolverDnDpotrf */
    gpu_p->d_devSync[gpuid] 	 = gpu_p->gpuPtr[gpuid]; 
    gpu_p->gpuPtr[gpuid] 		+= 2*sizeof(int)*gb_p->maxbatch; 
  
    /* cuSolver Cholesky initialization (get workspace size) */
    cusolverErr = cusolverDnDpotrf_bufferSize(Common->cusolverHandle[gpuid], CUBLAS_FILL_MODE_LOWER, gb_p->maxnscol, gpu_p->d_C[gpuid], gb_p->maxnsrow, &gb_p->work_size);
    if (cusolverErr != CUSOLVER_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU cusolverDnDpotrf_bufferSize failure");
    }

    /* Copy arrays to device from pinned memory */
    Ls = L->s;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;
    
    LsDim = gb_p->LsSize / sizeof(Int);
    ApDim = gb_p->ApSize / sizeof(Int);
    AiDim = gb_p->AiSize / sizeof(Int);
    AxDim = gb_p->AxSize / sizeof(double);
    
    /* Set pinned memory pointers */
    gpu_p->hostPtr[gpuid] 	 = Common->host_pinned_mempool[gpuid];

    /* Ax */
    h_Ax 				 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= gb_p->AxSize;
    
    /* Ls */
    h_Ls 				 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= gb_p->LsSize;

    /* Ap */
    h_Ap 				 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= gb_p->ApSize;   

    /* Ai */
    h_Ai 				 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= gb_p->AiSize;

    cudaStream_t copystream, memsetstream;
    cudaErr = cudaStreamCreate (&copystream);
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaStreamCreate(copystream)");
    cudaErr = cudaStreamCreate (&memsetstream);
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaStreamCreate(memsetstream)");
    
    /* Copy arrays to pinned memory from regular memory */  
#pragma omp parallel for num_threads(numThreads) private(i)
    for(i = 0; i < AxDim; i++) {
      h_Ax[i] = Ax[i];
    }
    cudaErr = cudaMemcpyAsync ( gpu_p->d_Ax[gpuid], h_Ax, gb_p->AxSize, cudaMemcpyHostToDevice, copystream);
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemcpy(d_Ax)");

#pragma omp parallel for num_threads(numThreads) private(i)
    for(i = 0; i < LsDim; i++){
      h_Ls[i] = Ls[i];
    }
    cudaErr = cudaMemcpyAsync ( gpu_p->d_Ls[gpuid], h_Ls, gb_p->LsSize, cudaMemcpyHostToDevice, copystream);
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemcpy(d_Ls)");	        

#pragma omp parallel for num_threads(numThreads) private(i)
    for(i = 0; i < ApDim; i++){
      h_Ap[i] = Ap[i];
    }
    cudaErr = cudaMemcpyAsync ( gpu_p->d_Ap[gpuid], h_Ap, gb_p->ApSize, cudaMemcpyHostToDevice, copystream);
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemcpy(d_Ap)");                

#pragma omp parallel for num_threads(numThreads) private(i)
    for(i = 0; i < AiDim; i++){
      h_Ai[i] = Ai[i]; 
    }
    cudaErr = cudaMemcpyAsync ( gpu_p->d_Ai[gpuid], h_Ai, gb_p->AiSize, cudaMemcpyHostToDevice, copystream);
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemcpy(d_Ai)");                
    
    /* Copy arrays to device from pinned memory */
    cudaErr = cudaMemsetAsync ( gpu_p->d_Lx[gpuid], 0, gb_p->LxSize, memsetstream);
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemset(d_Lx)");

    /* set pinned memory pointers */
    gpu_p->hostPtr[gpuid] 	 = Common->host_pinned_mempool[gpuid];

    /* type (double *) */

    gpu_p->h_ptrSuper[gpuid] 	 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= 3*gb_p->ptrSuperSize;

    gpu_p->h_ptrDesc[gpuid] 	 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= 6*gb_p->ptrDescSize;
  
    /* type (double) */

    gpu_p->h_Lx[gpuid] 		 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= gb_p->LxSize;
    
    /* type (int) */
    
    gpu_p->h_dimSuper[gpuid] 	 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= 13*gb_p->dimSuperSize;
    
    gpu_p->h_dimDesc[gpuid] 	 = gpu_p->hostPtr[gpuid];
    gpu_p->hostPtr[gpuid] 	+= 14*gb_p->dimDescSize;
    
    /* get device pointer to pinted buffer */
    cudaErr = cudaHostGetDevicePointer ( (void**)&gpu_p->h_pLx[gpuid], gpu_p->h_Lx[gpuid], 0);
    if (cudaErr) {
      printf("HostGetDevicePointer - cudaError:%s\n",cudaGetErrorString(cudaErr));
      ERROR ( CHOLMOD_GPU_PROBLEM, "GPU cudaHostGetDevicePointer");
    }
    
#ifdef USE_NVTX
    nvtxRangeEnd(range1);
#endif
  }

  return (1);  /* initialization successfull, useGPU = 1 */
  
}










/*
 *  Function:
 *    gpu_initialize_supernode_batch
 *
 *  Description:
 *    Initializes a batch of supernodes.
 *    Performs three tasks:
 *      1. sets device pointers
 *      2. creates map on device
 *      3. initializes Lx (factor) on device 
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_initialize_supernode_batch))
  (
   cholmod_common *Common,
   cholmod_global_pointers *gb_p,
   cholmod_gpu_pointers *gpu_p,         /* device pointers */
   Int n,
   Int maxsnsrow,			/* max nsrow for batch of supernodes */
   Int maxkdif,				/* max kdif for batch of supernodes */
   Int nzmax,				/* max nzmax for batch of supernodes */
   int strideSize,			/* stride for descendant lists */
   int nbatch,				/* size of batch */
   int gpuid			        /* gpu id */
   )
{

  /* local variables */
  int stream, maxbatch;
  struct cholmod_super_ptrs_t *h_super, *d_super;
  cudaError_t cudaStat;


  maxbatch = gb_p->maxbatch;
    
  /* 
   * Set pointers for batching arrays:    
   * sets pointers for storing dimensions and matrix addresses of 
   * supernodes and descendants, in the order they are batched.
   */

  /* set pointers for storing dimensions of batched supernodes */
  gpu_p->d_potrf[gpuid].n     	= gpu_p->d_dimSuper[gpuid] + 0*maxbatch;
  gpu_p->d_potrf[gpuid].lda   	= gpu_p->d_dimSuper[gpuid] + 1*maxbatch;
  gpu_p->d_trsm[gpuid].m   	= gpu_p->d_dimSuper[gpuid] + 2*maxbatch;
  gpu_p->d_trsm[gpuid].n   	= gpu_p->d_dimSuper[gpuid] + 3*maxbatch;
  gpu_p->d_trsm[gpuid].lda   	= gpu_p->d_dimSuper[gpuid] + 4*maxbatch;
  gpu_p->d_trsm[gpuid].ldb  	= gpu_p->d_dimSuper[gpuid] + 5*maxbatch;
  gpu_p->d_super[gpuid].s   	= gpu_p->d_dimSuper[gpuid] + 6*maxbatch;
  gpu_p->d_super[gpuid].k1   	= gpu_p->d_dimSuper[gpuid] + 7*maxbatch;
  gpu_p->d_super[gpuid].k2   	= gpu_p->d_dimSuper[gpuid] + 8*maxbatch;
  gpu_p->d_super[gpuid].psi  	= gpu_p->d_dimSuper[gpuid] + 9*maxbatch;
  gpu_p->d_super[gpuid].psx   	= gpu_p->d_dimSuper[gpuid] + 10*maxbatch;
  gpu_p->d_super[gpuid].nscol   = gpu_p->d_dimSuper[gpuid] + 11*maxbatch;
  gpu_p->d_super[gpuid].nsrow   = gpu_p->d_dimSuper[gpuid] + 12*maxbatch;

  /* set pointers for storing pointers of matrices of batched supernodes */
  gpu_p->d_potrf[gpuid].A 	= gpu_p->d_ptrSuper[gpuid] + 0*maxbatch;
  gpu_p->d_trsm[gpuid].A 	= gpu_p->d_ptrSuper[gpuid] + 1*maxbatch;
  gpu_p->d_trsm[gpuid].B 	= gpu_p->d_ptrSuper[gpuid] + 2*maxbatch;

  /* set pointers for storing dimensions of batched descendants */
  gpu_p->d_syrk[gpuid].n  	= gpu_p->d_dimDesc[gpuid] + 0*strideSize;
  gpu_p->d_syrk[gpuid].k  	= gpu_p->d_dimDesc[gpuid] + 1*strideSize;
  gpu_p->d_syrk[gpuid].lda  	= gpu_p->d_dimDesc[gpuid] + 2*strideSize;
  gpu_p->d_syrk[gpuid].ldc  	= gpu_p->d_dimDesc[gpuid] + 3*strideSize;
  gpu_p->d_gemm[gpuid].m  	= gpu_p->d_dimDesc[gpuid] + 4*strideSize;
  gpu_p->d_gemm[gpuid].n  	= gpu_p->d_dimDesc[gpuid] + 5*strideSize;
  gpu_p->d_gemm[gpuid].k  	= gpu_p->d_dimDesc[gpuid] + 6*strideSize;
  gpu_p->d_gemm[gpuid].lda  	= gpu_p->d_dimDesc[gpuid] + 7*strideSize;
  gpu_p->d_gemm[gpuid].ldb  	= gpu_p->d_dimDesc[gpuid] + 8*strideSize;
  gpu_p->d_gemm[gpuid].ldc  	= gpu_p->d_dimDesc[gpuid] + 9*strideSize;
  gpu_p->d_desc[gpuid].ndrow1  	= gpu_p->d_dimDesc[gpuid] + 10*strideSize;
  gpu_p->d_desc[gpuid].ndrow2  	= gpu_p->d_dimDesc[gpuid] + 11*strideSize;
  gpu_p->d_desc[gpuid].pdi1  	= gpu_p->d_dimDesc[gpuid] + 12*strideSize;
  gpu_p->d_desc[gpuid].s  	= gpu_p->d_dimDesc[gpuid] + 13*strideSize;

  /* set pointers for storing pointers of matrices of batched descendants */
  gpu_p->d_syrk[gpuid].A 	= gpu_p->d_ptrDesc[gpuid] + 0*strideSize;
  gpu_p->d_syrk[gpuid].C 	= gpu_p->d_ptrDesc[gpuid] + 1*strideSize;
  gpu_p->d_gemm[gpuid].A 	= gpu_p->d_ptrDesc[gpuid] + 2*strideSize;
  gpu_p->d_gemm[gpuid].B 	= gpu_p->d_ptrDesc[gpuid] + 3*strideSize;
  gpu_p->d_gemm[gpuid].C 	= gpu_p->d_ptrDesc[gpuid] + 4*strideSize;
  gpu_p->d_desc[gpuid].C 	= gpu_p->d_ptrDesc[gpuid] + 5*strideSize;




  /* set dimensions*/
  h_super    = &gpu_p->h_super[gpuid];
  d_super    = &gpu_p->d_super[gpuid];




  /* create map for batch of supernodes */
  createMapOnDevice_batch ( gpu_p->d_Map[gpuid],
                	    gpu_p->d_Ls[gpuid],
                      	    d_super->psi,
                       	    d_super->nsrow,
                      	    maxsnsrow,
			    n,
                      	    nbatch,
                      	    &(Common->gpuStream[gpuid][0]));

  if (cudaGetLastError()!=cudaSuccess) {
    printf("error: %s\n",cudaGetErrorString(cudaGetLastError()));
    ERROR (CHOLMOD_GPU_PROBLEM, "GPU createMapOnDevice_batch\n") ;
  }




  /* initialize Lx (factor) for batch of supernoeds */
  initLxonDevice_batch( gpu_p->d_Lx[gpuid],
                        gpu_p->d_Ax[gpuid],
                        gpu_p->d_Ap[gpuid],
                        gpu_p->d_Ai[gpuid],
                        gpu_p->d_Map[gpuid],
                        d_super->nsrow,
		        d_super->psx,
		        d_super->k1,
		        d_super->k2,
                        nzmax,
                        maxkdif,
                        n,
		        nbatch,
                        &(Common->gpuStream[gpuid][0]));

  if (cudaGetLastError()!=cudaSuccess) {
    printf("error: %s\n",cudaGetErrorString(cudaGetLastError()));
    ERROR (CHOLMOD_GPU_PROBLEM, "GPU initLxonDevice_batch\n") ;
  }




  /* synchronize stream */
  cudaStat = cudaStreamSynchronize (Common->gpuStream[gpuid][0]) ;
  if(cudaStat) {
    ERROR (CHOLMOD_GPU_PROBLEM, "GPU gpu_initialize_supernode_batch\n") ;
  }
 
}










/*
 *  Function:
 *    gpu_updateC_batch
 *
 *  Description:
 *    Computes the schur compliment (DSYRK & DGEMM) of a batch of supernodes and maps
 *    it back (addUpdate). 
 *    Performs three tasks:
 *      1. DSYRK
 *        a. cuBlas calls (n > D_N || k > D_K)
 *        b. batched call (otherwise)
 *
 *      2. DGEMM
 *        a. cuBlas calls (m > D_M || n > D_N || k > D_K)
 *        b. batched call (otherwise)
 *
 *      3. addUpdate (map schur comp. to supernode)  
 *        a. update one supernode (ndrow1 > D_ROW1 || ndrow2 > D_ROW2)
 *        b. batched call (otherwise)
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_updateC_batch))
  (
   cholmod_common *Common,
   cholmod_global_pointers *gb_p,
   cholmod_gpu_pointers *gpu_p,
   cholmod_cpu_pointers *cpu_p,
   cholmod_tree_pointers *tree_p,
   cholmod_profile_pointers *prof_p,
   int n,   
   int subtree,
   int level,
   int nbatch,			/* size of batch (# descendants) */
   int syrk_count,		/* # syrk calls sent to cuBlas */
   int gemm_count,		/* # gemm calls sent to cuBlas */
   int update_count,		/* # addUpdate calls to stream seperately */
   int gpuid,
   int flag,
   Int numSuper,
   Int start,
   Int end,
   int *max_dim,            /* max m,n dimensions for batch */
   Int *LpxSub,
   Int *Lpx
   )
{
  /* local variables */
  int i;
  int *ndrow1, *ndrow2, *ndrow3, *ndrow, *ndcol, *pdi1, *list, *nsrow, *psx;
  double alpha, beta, tstart1;
  double **Aptr, **Bptr, **Cptr, **C2ptr, *tend, *gemm_time, *syrk_time;
  struct cholmod_syrk_ptrs_t *h_syrk, *d_syrk;
  struct cholmod_gemm_ptrs_t *h_gemm, *d_gemm;
  struct cholmod_desc_ptrs_t *h_desc, *d_desc;
  struct cholmod_super_ptrs_t *h_super, *d_super;
  cublasStatus_t cublasStatus ;
  cudaError_t cudaStat;
  
  alpha  = -1.0 ;
  beta   = 0.0 ;





  /* set profile pointers */
  tend          = prof_p->f_end[gpuid];
  syrk_time     = prof_p->syrk_time[gpuid];
  gemm_time     = prof_p->gemm_time[gpuid];





  /* 
   * Perform DSYRK  
   * Compute dsyrk for all descendants 
   * of a batch of supernodes.
   */
  TIMER_START1(tstart1);

  /* set dimensions */
  h_syrk    = &gpu_p->h_syrk[gpuid];
  d_syrk    = &gpu_p->d_syrk[gpuid];

  /* loop over 'large' syrk's */
  for(i = 0; i < syrk_count; i++) {

    /* set cublas stream */
    cublasStatus = cublasSetStream (Common->cublasHandle[gpuid], Common->gpuStream[gpuid][i%CHOLMOD_DEVICE_STREAMS]) ;
    if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU cublasSetStream failure");
    }

    /* dsyrk on cuBlas */
    cublasDsyrk ( Common->cublasHandle[gpuid],
                  CUBLAS_FILL_MODE_LOWER,
                  CUBLAS_OP_N,
                  h_syrk->n[i],
                  h_syrk->k[i],
                  &alpha,
                  h_syrk->A[i],
                  h_syrk->lda[i],
                  &beta,
                  h_syrk->C[i],
                  h_syrk->ldc[i]);
  }


  /* check if any syrk's left for batching */
  if( (nbatch - syrk_count) > 0 ) {

    /* dsyrk on batched kernels */   
    dsyrk_custom_simple_1block_batch( Common->gpuStream[gpuid][0],
                                      CUBLAS_FILL_MODE_LOWER,
                                      CUBLAS_OP_N,
                                      &d_syrk->n[syrk_count],
                                      &d_syrk->k[syrk_count],
                                      &alpha,
                                      &d_syrk->A[syrk_count],
                                      &d_syrk->lda[syrk_count],
                                      &beta,
                                      &d_syrk->C[syrk_count],
                                      &d_syrk->ldc[syrk_count],
                                      nbatch-syrk_count);
  }


  /* synchronize streams */
  for(i = 0; i < CHOLMOD_DEVICE_STREAMS; i++){
    cudaStat = cudaStreamSynchronize (Common->gpuStream[gpuid][i]) ;
    if (cudaStat) { 
      printf("error: %s\n",cudaGetErrorString(cudaGetLastError()));
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU dsyrk_custom_simple_1block_batch") ; 
    }
  }
  TIMER_END1(tstart1,syrk_time,0);
  TIMER_END1(tstart1,syrk_time,1);



  /* 
   * Perform DGEMM 
   * Compute dgemm for all descendants 
   * of a batch of supernodes.
   */
  TIMER_START1(tstart1);

  /* set dimensions */
  h_gemm  = &gpu_p->h_gemm[gpuid];
  d_gemm  = &gpu_p->d_gemm[gpuid];


  /* loop over 'large' gemm's */
  for(i = 0; i < gemm_count; i++) {

    /* set cublas stream */
    cublasStatus = cublasSetStream (Common->cublasHandle[gpuid], Common->gpuStream[gpuid][i%CHOLMOD_DEVICE_STREAMS]) ;
    if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU cublasSetStream failure");
    }

    /* dgemm on cuBlas */
    cublasDgemm ( Common->cublasHandle[gpuid],
                  CUBLAS_OP_N, CUBLAS_OP_T,
                  h_gemm->m[i],
                  h_gemm->n[i],
                  h_gemm->k[i],
                  &alpha,
                  h_gemm->A[i],
                  h_gemm->lda[i],
                  h_gemm->B[i],
                  h_gemm->ldb[i],
                  &beta,
                  h_gemm->C[i],
                  h_gemm->ldc[i]);
  }

  
  /* check if any gemm's left for batching */
  if( (nbatch - gemm_count) > 0 ) {   

    /* dgemm on batched kernels */
    dgemm_custom_simple_1block_batch( Common->gpuStream[gpuid][0],
                                      CUBLAS_OP_N, CUBLAS_OP_T,
                                      &d_gemm->m[gemm_count],
                                      &d_gemm->n[gemm_count],
                                      &d_gemm->k[gemm_count],
                                      &alpha,
                                      &d_gemm->A[gemm_count],
                                      &d_gemm->lda[gemm_count],
                                      &d_gemm->B[gemm_count],
                                      &d_gemm->ldb[gemm_count],
                                      &beta,
                                      &d_gemm->C[gemm_count],
                                      &d_gemm->ldc[gemm_count],
                                      nbatch-gemm_count);
  } 


  /* copy supernode from pinned to regular memory - only at last level */
  if ( level == tree_p->supernode_num_levels[subtree]-1 ) {

  TEMPLATE2 (CHOLMOD(gpu_copy_supernode))( Common,
                                           gpu_p,
					   cpu_p,
					   tree_p,
                                           subtree,
                                           level,
                                           gpuid,
                                           1,
                                           numSuper,
                                           start,
                                           end,
                                           LpxSub,
                                           Lpx);

  }


  /* synchronize streams */
  for(i = 0; i < CHOLMOD_DEVICE_STREAMS; i++){
    cudaStat = cudaStreamSynchronize (Common->gpuStream[gpuid][i]) ;
    if (cudaStat) { 
      printf("error: %s\n",cudaGetErrorString(cudaGetLastError()));
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU dgemm_custom_simple_1block_batch") ; 
    }
  }
  TIMER_END1(tstart1,gemm_time,0);
  TIMER_END1(tstart1,gemm_time,1);





  /* 
   * Perform addUpdate 
   * Add schur compliment for all descendants
   * of a batch of supernodes.  
   */
  TIMER_START1(tstart1);

  /* set dimensions */
  h_desc  = &gpu_p->h_desc[gpuid];
  d_desc  = &gpu_p->d_desc[gpuid];
  h_super  = &gpu_p->h_super[gpuid];
  d_super  = &gpu_p->d_super[gpuid];



  /* loop over 'large' addUpdate's */
  for(i = 0; i < update_count; i++) {

    /* mapping (to factor Lx) for each descendants */
    addUpdateOnDevice_large( gpu_p->d_Lx[gpuid],
                       	     &d_desc->C[i],
                       	     gpu_p->d_Map[gpuid],
                             gpu_p->d_Ls[gpuid],
                             h_desc->pdi1[i],
                             h_desc->ndrow1[i],
                             h_desc->ndrow2[i],
                             h_super->psx[h_desc->s[i]],
                             h_super->nsrow[h_desc->s[i]],
		             (int)((n+1)*h_desc->s[i]),
                             &(Common->gpuStream[gpuid][i%CHOLMOD_DEVICE_STREAMS]));
  }

  /* check if any addUpdate's left for batching */
  if( (nbatch - update_count) > 0 ) {

    /* mapping (to factor Lx) for descendants in a batch of supernodes */
    addUpdateOnDevice_batch( gpu_p->d_Lx[gpuid],
                             &d_desc->C[update_count],
                             gpu_p->d_Map[gpuid],
                             gpu_p->d_Ls[gpuid],
                             n,
                             &d_desc->pdi1[update_count],
                             &d_desc->ndrow1[update_count],
                             &d_desc->ndrow2[update_count],
                             d_super->psx,
                             d_super->nsrow,
                             &d_desc->s[update_count],
                             max_dim,
                             nbatch-update_count,
                             &(Common->gpuStream[gpuid][0]));

  }


  /* synchronize streams */
  for(i = 0; i < CHOLMOD_DEVICE_STREAMS; i++){
    cudaStat = cudaStreamSynchronize (Common->gpuStream[gpuid][i]) ;
    if (cudaStat) {
      printf("error: %s\n",cudaGetErrorString(cudaGetLastError()));
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU addUpdateOnDevice_batch") ;
    }
  }
  TIMER_END1(tstart1,tend,4);

}










/*
 *  Function:
 *    gpu_lower_potrf_batch
 *
 *  Description:
 *    computes cholesky factoriation for a batch of supernodes. 
 *    Performs one task:
 *      1. DPOTRF
 *        a. cuSolver calls (n < S_N)
 *        b. batched call
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_lower_potrf_batch))
  (
   cholmod_common *Common,
   cholmod_global_pointers *gb_p,
   cholmod_gpu_pointers *gpu_p,		/* device pointers */
   cholmod_profile_pointers *prof_p,
   int nbatch,				/* batch size (# supernodes) */
   int potrf_count,			/* # potrf calls sent to cuSolver */
   int gpuid				/* gpu id */
   )
{
  /* local variables */
  int i, info;
  int *nscol2, *nsrow;
  double tstart1;
  double **Aptr, *tend, *potrf_time;
  struct cholmod_potrf_ptrs_t *h_potrf, *d_potrf;
  cudaError_t cudaStat;
  cusolverStatus_t cusolverErr;
    

    

  /* set profile pointers */
  tend          = prof_p->f_end[gpuid];
  potrf_time     = prof_p->potrf_time[gpuid];




  /* 
   * Perform batched DPOTRF 
   * Compute potrf for all supernodes
   * in a batch  
   */
  TIMER_START1(tstart1);
   
  /* set dimensions */
  h_potrf  = &gpu_p->h_potrf[gpuid];
  d_potrf  = &gpu_p->d_potrf[gpuid];

  /* loop over 'large' potrf's */
  for(i=0; i<potrf_count; i++) {


    /* set cuSolver stream */
    cusolverErr = cusolverDnSetStream (Common->cusolverHandle[gpuid], Common->gpuStream[gpuid][i%CHOLMOD_DEVICE_STREAMS]) ;
    if (cusolverErr != CUSOLVER_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU cusolverDnSetStream failure");
    }
  
    /* potrf on cuSolver */
    cusolverDnDpotrf( Common->cusolverHandle[gpuid], 
		      CUBLAS_FILL_MODE_LOWER, 
		      h_potrf->n[i],
		      h_potrf->A[i],
		      h_potrf->lda[i],
		      &gpu_p->d_devSync[gpuid][2*i], 
		      gb_p->work_size, 
		      gpu_p->d_info[gpuid]);

  }


  /* check if any potrf's left for batching */
  if( (nbatch - potrf_count) > 0 ) {
  
    /* potrf on custom batched kernels */
    dpotrf_custom_simple_1block_batch ( Common->gpuStream[gpuid][0],
                                        CUBLAS_FILL_MODE_LOWER,
                                        &d_potrf->n[i],
                                        &d_potrf->A[i],
                                        &d_potrf->lda[i],
                                        &info,
                                        nbatch-potrf_count) ;

  } 


  /* synchronize streams */
  for(i = 0; i < CHOLMOD_DEVICE_STREAMS; i++){
    cudaStat = cudaStreamSynchronize (Common->gpuStream[gpuid][i]) ;
    if (cudaStat) { 
      printf("error: %s\n",cudaGetErrorString(cudaGetLastError()));
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU dpotrf_custom_simple_1block_batch") ; 
    }
  }
  TIMER_END1(tstart1,potrf_time,0);
  TIMER_END1(tstart1,potrf_time,1);

}










/*
 *  Function:
 *    gpu_triangular_solve_batch
 *
 *  Description:
 *    Computes triangular solve for a batch of supernodes.
 *    Performs one task:
 *      1. DTRSM
 *        a. cuBlas calls (m < S_M || n < S_N)
 *        b. batched call (otherwise)
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_triangular_solve_batch))
  (
   cholmod_common *Common,
   cholmod_global_pointers *gb_p,
   cholmod_gpu_pointers *gpu_p,		/* device pointers */
   cholmod_profile_pointers *prof_p,
   int nbatch,				/* batch size (# supernodes) */
   int trsm_count,			/* # trsm calls sent to cuBlas */
   int gpuid				/* gpu id */
   )
{
  /* local variables */
  int i;
  int *nsrow2, *nscol2, *nsrow, *psx;
  double alpha, tstart1;
  double **Aptr, **Bptr, *tend, *trsm_time;
  struct cholmod_trsm_ptrs_t *h_trsm, *d_trsm;
  cublasStatus_t cublasStatus ;
  cudaError_t cudaStat ;
   
  alpha  = 1.0 ;




  /* set profile pointers */
  tend          = prof_p->f_end[gpuid];
  trsm_time     = prof_p->trsm_time[gpuid];




  /* 
   * Perform batched DTRSM 
   * Compute trsm for all supernodes
   * in a batch
   */
  TIMER_START1(tstart1);

  /* set dimensions */
  h_trsm  = &gpu_p->h_trsm[gpuid];
  d_trsm  = &gpu_p->d_trsm[gpuid];


  /* loop over 'large' trsm's */
  for(i=0; i<trsm_count; i++) {

    /* set cublas stream */
    cublasStatus = cublasSetStream (Common->cublasHandle[gpuid], Common->gpuStream[gpuid][i%CHOLMOD_DEVICE_STREAMS]) ;
    if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU cublasSetStream failure");
    }

    /* trsm on cuBlas */
    cublasDtrsm ( Common->cublasHandle[gpuid],
                  CUBLAS_SIDE_RIGHT,
                  CUBLAS_FILL_MODE_LOWER,
                  CUBLAS_OP_T,
                  CUBLAS_DIAG_NON_UNIT,
                  h_trsm->m[i],
                  h_trsm->n[i],
                  &alpha,
                  h_trsm->A[i],
                  h_trsm->lda[i],
                  h_trsm->B[i],
                  h_trsm->ldb[i]);
  }


  /* check if any trsm's left for batching */
  if( (nbatch - trsm_count) > 0 ) {

    /* trsm on custom batched kernels */
    dtrsm_custom_simple_1block_batch ( Common->gpuStream[gpuid][0],
                                       CUBLAS_SIDE_RIGHT,
                                       CUBLAS_FILL_MODE_LOWER,
                                       CUBLAS_OP_T,
                                       CUBLAS_DIAG_NON_UNIT,
                                       &d_trsm->m[i],
                                       &d_trsm->n[i],
                                       &alpha,
                                       &d_trsm->A[i],
                                       &d_trsm->lda[i],
                                       &d_trsm->B[i],
                                       &d_trsm->ldb[i],
                                       nbatch-trsm_count);
  } 


  /* synchronize streams */
  for(i = 0; i < CHOLMOD_DEVICE_STREAMS; i++){
    cudaStat = cudaStreamSynchronize (Common->gpuStream[gpuid][i]) ;
    if (cudaStat) {
      printf("error: %s\n",cudaGetErrorString(cudaGetLastError())); 
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU dtrsm_custom_simple_1block_batch") ; 
    }
  }
  TIMER_END1(tstart1,trsm_time,0);
  TIMER_END1(tstart1,trsm_time,1);

}










/*
 *  Function: 
 *    gpu_copy_supernode2
 *
 *  Description:
 *    Copies batch of supernodes from device to pinned memory.
 *    Performs zero-copies.
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_copy_supernode2))
  (
   cholmod_common *Common,
   cholmod_global_pointers *gb_p,
   cholmod_gpu_pointers *gpu_p,		/* device pointers */
   int nbatch,				/* batch size (# supernodes) */
   int maxnsrownscol,			/* max nsrow*nscol in batch */
   int gpuid				/* gpu id */
   )
{
  /* local variables */
  int *nscol2, *nsrow, *psx;
  struct cholmod_super_ptrs_t *h_super, *d_super;
  cudaError_t cudaStat ;


  /* set dimensions */
  h_super  = &gpu_p->h_super[gpuid];
  d_super  = &gpu_p->d_super[gpuid];


  /*
   * Perform copy factor
   * Copies a batch of supernodes from
   * device to pinned memory.
   */
  copyLx_small ( gpu_p->h_pLx[gpuid],
                 gpu_p->d_Lx[gpuid],
                 d_super->psx,
                 d_super->nsrow,
                 d_super->nscol,
                 nbatch,
                 maxnsrownscol,
                 &Common->gpuStream[gpuid][0]);

}










/*
 *  Function:
 *    gpu_copy_supernode
 *
 *  Description:
 *    copies batch of supernodes from pinned to regular memory 
 *    Performs two tasks:
 *      1. copy supernodes up to last level (flag=1)
 *      2. copy supernodes in last level (flag=2)
 *
 *  Note: 
 *    in second case (flag=2) better to set omp threads across
 *    work over single supernode, as there are few (generally
 *    only one) supernodes in the final level. 
 *        
 */
void TEMPLATE2 (CHOLMOD (gpu_copy_supernode))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,		/* device pointers */
   cholmod_cpu_pointers *cpu_p,
   cholmod_tree_pointers *tree_p,
   int subtree,				/* subtree */
   int level,				/* level */
   int gpuid,				/* gpu id */
   int flag,				/* flag (1: copy up to last level, 2: copy last level) */
   Int numSuper,			/* # supernodes */
   Int start,				/* first supernode in level  */
   Int end,				/* last supernode in level */
   Int *LpxSub,				/* sub-factor in pinned memory */
   Int *Lpx				
   )
{
  /* local variables */
  Int i, j, s, nscol, nsrow, offset1, offset2;
  cudaError_t cudaStat ;
  int numThreads;



  numThreads = Common->ompNumThreads;




  /* case 1: copy all supernodes up to last level */
  if(flag==1) {

    /* split various supernodes amongst omp threads */
    #pragma omp parallel for num_threads(numThreads) schedule(static,1)
    for(i = 0; i < numSuper; i++) {
      s = tree_p->supernode_subtree[tree_p->supernode_subtree_ptrs[subtree] + i];
      nscol = cpu_p->Super [s+1] - cpu_p->Super [s] ;
      nsrow = cpu_p->Lpi[s+1] - cpu_p->Lpi[s] ;
      const int ssize = nscol*nsrow;
 
      offset1 = ((Int*)Lpx)[s];
      offset2 = LpxSub[s];
      double* pLx = &(cpu_p->Lx[offset1]);
      const double* phLx = &(gpu_p->h_Lx[gpuid][offset2]);

      /* copy supernode */
      for(j = 0; j < ssize; j++) {
        *pLx++ = *phLx++;
      }

    } /* end loop over supernodes */
  } /* end case 1*/




  /* case 2: copy all supernodes in last level */
  if(flag==2) {

    /* loop over supernodes in last level */
    for(i = 0; i < (end - start); i++) {

      s = tree_p->supernode_levels[start + i];
      nscol = cpu_p->Super [s+1] - cpu_p->Super [s] ;
      nsrow = cpu_p->Lpi[s+1] - cpu_p->Lpi[s] ;
      const int ssize = nscol*nsrow;

      offset1 = ((Int*)Lpx)[s];
      offset2 = LpxSub[s];
      double* pLx = &(cpu_p->Lx[offset1]);
      const double* phLx = &(gpu_p->h_Lx[gpuid][offset2]);

      /* copy supernode - split single supernode amongst omp threads  */
      #pragma omp parallel for private(j)
      for(j = 0; j < ssize; j++) {
        pLx[j] = phLx[j];
      }

    } /* end loop over supernodes */
  }
   
}




#endif


/*
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
*/












