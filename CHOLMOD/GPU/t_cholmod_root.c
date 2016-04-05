/* ========================================================================== */
/* === GPU/t_cholmod_root =================================================== */
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
 *   t_cholmod_root
 *
 * Description:
 *   Contains functions for root algorithm.
 *
 */


/* includes */
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <string.h>
#include <time.h>


/* macros */
#undef L_ENTRY
#ifdef REAL
#define L_ENTRY 1
#else
#define L_ENTRY 2
#endif










/*
 * Function:
 *   gpu_init_root
 *
 * Description:
 *   Performs required initialization for GPU computing.
 *   Returns 0 if there is an error, disabling GPU computing (useGPU = 0)
 *
 */
int TEMPLATE2 (CHOLMOD (gpu_init_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,
   cholmod_factor *L,
   Int *Lpi,
   Int nsuper,
   Int n,
   int gpuid
   )
{
  /* local variables */
  Int i, k, nls;
  cublasStatus_t cublasError;
  cudaError_t cudaErr;


#ifdef USE_NVTX
  nvtxRangeId_t range1 = nvtxRangeStartA("gpu_init_root");
#endif


  /* set cuda device */
  cudaSetDevice(gpuid);


  /* compute nls */
  nls =  Lpi[nsuper]-Lpi[0];


  /* make buffer size is large enough */
  if ( (nls+2*n+4)*sizeof(Int) > Common->devBuffSize ) {
    ERROR (CHOLMOD_GPU_PROBLEM,"\n\n"
           "GPU Memory allocation error.  Ls, Map and RelativeMap exceed\n"
           "devBuffSize.  It is not clear if this is due to insufficient\n"
           "device or host memory or both.  You can try:\n"
           "     1) increasing the amount of GPU memory requested\n"
           "     2) reducing CHOLMOD_NUM_HOST_BUFFERS\n"
           "     3) using a GPU & host with more memory\n"
           "This issue is a known limitation and should be fixed in a \n"
           "future release of CHOLMOD.\n") ;
    return (0) ;
  }




  /* set gpu memory pointers */

  /* type double */
  gpu_p->d_Lx_root[gpuid][0] = Common->dev_mempool[gpuid];
  gpu_p->d_Lx_root[gpuid][1] = Common->dev_mempool[gpuid] + Common->devBuffSize;
  gpu_p->d_C_root[gpuid] = Common->dev_mempool[gpuid] + 2*Common->devBuffSize;
  gpu_p->d_A_root[gpuid][0] = Common->dev_mempool[gpuid] + 3*Common->devBuffSize;
  gpu_p->d_A_root[gpuid][1] = Common->dev_mempool[gpuid] + 4*Common->devBuffSize;

  /* type Int */
  gpu_p->d_Ls_root[gpuid] = Common->dev_mempool[gpuid] + 5*Common->devBuffSize;
  gpu_p->d_Map_root[gpuid] = Common->dev_mempool[gpuid] + 5*Common->devBuffSize + (nls+1)*sizeof(Int);
  gpu_p->d_RelativeMap_root[gpuid] = Common->dev_mempool[gpuid] + 5*Common->devBuffSize + (nls+1)*sizeof(Int) + (n+1)*sizeof(Int);




  /* copy Ls and Lpi to device */
  cudaErr = cudaMemcpy ( gpu_p->d_Ls_root[gpuid], L->s, nls*sizeof(Int), cudaMemcpyHostToDevice );
  CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemcpy(d_Ls_root)");




  /* set pinned memory pointers */
  gpu_p->h_Lx_root[gpuid][0] = (double*)(Common->host_pinned_mempool[gpuid]);

  for (k = 1; k < CHOLMOD_HOST_SUPERNODE_BUFFERS; k++) {
    gpu_p->h_Lx_root[gpuid][k] = (double*)((char *)(Common->host_pinned_mempool[gpuid]) + k*Common->devBuffSize);
  }


#ifdef USE_NVTX
  nvtxRangeEnd(range1);
#endif


  return (1);  /* initialization successfull, useGPU = 1 */

}










/*
 * Function:
 *   gpu_reorder_descendants_root
 *
 * Description:
 *   Reorders descendants in a supernode by size (ndrow2*ndcol) 
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_reorder_descendants_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,
   Int *Lpi,
   Int *Lpos,
   Int *Super,
   Int *Head,
   Int *tail,
   Int *Next,
   Int *Previous,
   Int *ndescendants,
   Int *mapCreatedOnGpu,
   Int locals,
   int gpuid
   )
{
  /* local variables */
  Int d, k, p, kd1, kd2, ndcol, ndrow2, pdi, pdend, pdi1, nextd, dnext, n_descendant = 0;
  int previousd, nreverse = 1, numThreads;
  double score;

  /* store GPU-eligible descendants in h_Lx[0] */
  struct cholmod_descendant_score_t* scores = (struct cholmod_descendant_score_t*) gpu_p->h_Lx_root[gpuid][0];



  /* initialize variables */
  d = Head[locals];
  *mapCreatedOnGpu = 0;
  numThreads	= Common->ompNumThreads;



  /* loop until reach last descendant in supernode */
  while ( d != EMPTY ) {

      /* get dimensions for the current descendant */
      kd1 = Super [d] ;       /* d contains cols kd1 to kd2-1 of L */
      kd2 = Super [d+1] ;
      ndcol = kd2 - kd1 ;     /* # of columns in all of d */
      pdi = Lpi [d] ;         /* pointer to first row of d in Ls */
      pdend = Lpi [d+1] ;     /* pointer just past last row of d in Ls */
      p = Lpos [d] ;          /* offset of 1st row of d affecting s */
      pdi1 = pdi + p ;        /* ptr to 1st row of d affecting s in Ls */
      ndrow2 = pdend - pdi1;
     
      nextd = Next[d];

      /* compute the descendant's rough flops 'score' */
      score = ndrow2 * ndcol;
      if ( (ndrow2*L_ENTRY >= CHOLMOD_ND_ROW_LIMIT) && (ndcol*L_ENTRY >= CHOLMOD_ND_COL_LIMIT) ) {
        score += Common->devBuffSize;
      }

      /* store descendant in list */
      scores[n_descendant].score = score;
      scores[n_descendant].d = d;
      n_descendant++;
 
      d = nextd;

  } 



  /* sort the GPU-eligible supernodes in descending size (flops) */
  qsort ( scores, n_descendant, sizeof(struct cholmod_descendant_score_t), (__compar_fn_t) CHOLMOD(score_comp) );



  /* place sorted data back in descendant supernode linked list */
  if ( n_descendant > 0 ) {

    Head[locals] = scores[0].d;
    if ( n_descendant > 1 ) {

      #pragma omp parallel for num_threads(numThreads) if (n_descendant > 64)
      for ( k=1; k<n_descendant; k++ ) {
        Next[scores[k-1].d] = scores[k].d;
      }

    }

    Next[scores[n_descendant-1].d] = -1;
  }



  /* reverse the first CHOLMOD_HOST_SUPERNODE_BUFFERS descendants to better hide PCIe communications */
  if ( Head[locals] != EMPTY && Next[Head[locals]] != EMPTY ) {

    previousd = Head[locals];
    d = Next[Head[locals]];

    /* loop through the first CHOLMOD_HOST_SUPERNODE_BUFFERS descendants */
    while ( d!=EMPTY && nreverse < CHOLMOD_HOST_SUPERNODE_BUFFERS ) {

      /* get descendant dimensions */
      kd1 = Super [d] ;       /* d contains cols kd1 to kd2-1 of L */
      kd2 = Super [d+1] ;
      ndcol = kd2 - kd1 ;     /* # of columns in all of d */
      pdi = Lpi [d] ;         /* pointer to first row of d in Ls */
      pdend = Lpi [d+1] ;     /* pointer just past last row of d in Ls */
      p = Lpos [d] ;          /* offset of 1st row of d affecting s */
      pdi1 = pdi + p ;        /* ptr to 1st row of d affecting s in Ls */
      ndrow2 = pdend - pdi1;

      nextd = Next[d];

      nreverse++;

      /* place descendant at the front of the list */
      if ( (ndrow2*L_ENTRY >= CHOLMOD_ND_ROW_LIMIT) && (ndcol*L_ENTRY >= CHOLMOD_ND_COL_LIMIT) ) {
        Next[previousd] = Next[d];
        Next[d] = Head[locals];
        Head[locals] = d;
      }
      else {
        previousd = d;
      }

      d = nextd;
    } /* end while loop */

  } 



  /* create a 'previous' list so we can traverse backwards */
  n_descendant = 0;

  if ( Head[locals] != EMPTY ) {

    Previous[Head[locals]] = EMPTY;

    /* loop over descendants */
    for (d = Head [locals] ; d != EMPTY ; d = dnext) {

      n_descendant++;
      dnext = Next[d];

      if ( dnext != EMPTY ) {
        Previous[dnext] = d;
      }
      else {
        *tail = d;
      }

    } /* end loop over descendants */

  }


  /* store descendant dimension */
  *ndescendants = n_descendant;

}










/*
 *  Function:
 *    gpu_initialize_supernode_root
 *
 *  Description:
 *    Initializes a supernode.
 *    Performs two tasks:
 *      1. clears A buffer for assembly
 *      2. creates map on device
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_initialize_supernode_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,
   Int nscol,
   Int nsrow,
   Int psi,
   int gpuid
   )
{
  /* local variables */
  cudaError_t cudaErr;


  /* initialize the device supernode assemby memory to zero */
  cudaErr = cudaMemset ( gpu_p->d_A_root[gpuid][0], 0, nscol*nsrow*L_ENTRY*sizeof(double) );
  CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemset(d_A_root)");


  /* create the map for supernode on the device */
  createMapOnDevice ( (Int *)(gpu_p->d_Map_root[gpuid]), (Int *)(gpu_p->d_Ls_root[gpuid]), psi, nsrow );
  cudaErr = cudaGetLastError();
  CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"createMapOnDevice error");

}










/*
 *  Function:
 *    gpu_updateC_root
 *
 *  Description:
 *    Computes the schur compliment (DSYRK & DGEMM) of a batch of supernodes and maps
 *    it back (addUpdate).
 *    Performs three tasks:
 *      1. DSYRK
 *      2. DGEMM
 *      3. addUpdate (map schur comp. to supernode)
 */
int TEMPLATE2 (CHOLMOD (gpu_updateC_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,
   double *Lx,
   double *C,
   Int ndrow1,         
   Int ndrow2,
   Int ndrow,          
   Int ndcol,          
   Int nsrow,
   Int pdx1,           
   Int pdi1,
   int gpuid
   )
{
  /* local variables */
  int icol, irow, iHostBuff, iDevBuff, numThreads;
  Int ndrow3;
  double alpha, beta; 
  double *devPtrLx, *devPtrC;
  cublasStatus_t cublasStatus;
  cudaError_t cudaErr;


  numThreads = Common->ompNumThreads;




  /* early exit if descendant too small for cuBlas */
  if ( (ndrow2*L_ENTRY < CHOLMOD_ND_ROW_LIMIT) || (ndcol*L_ENTRY <  CHOLMOD_ND_COL_LIMIT) ) 
  {
    return (0) ;
  }



  /* initialize variables */
  ndrow3 = ndrow2 - ndrow1 ;
  alpha  = 1.0 ;
  beta   = 0.0 ;

  iHostBuff = (Common->ibuffer[gpuid])%CHOLMOD_HOST_SUPERNODE_BUFFERS;
  iDevBuff = (Common->ibuffer[gpuid])%2;

  /* initialize poitners */
  devPtrLx = (double *)(gpu_p->d_Lx_root[gpuid][iDevBuff]);
  devPtrC = (double *)(gpu_p->d_C_root[gpuid]);




  /*
   * Copy Lx to the device:
   * First copy to pinned buffer, then to the device for
   * better H2D bandwidth.
   */
  /* copy host data to pinned buffer */
#pragma omp parallel for num_threads(numThreads) if (ndcol > 32)
  for ( icol=0; icol<ndcol; icol++ ) {
    for ( irow=0; irow<ndrow2*L_ENTRY; irow++ ) {
      gpu_p->h_Lx_root[gpuid][iHostBuff][icol*ndrow2*L_ENTRY+irow] =
        Lx[pdx1*L_ENTRY+icol*ndrow*L_ENTRY + irow];
    }
  }



  /* copy pinned buffer to device */
  cudaErr = cudaMemcpyAsync ( devPtrLx,
                              gpu_p->h_Lx_root[gpuid][iHostBuff],
                              ndrow2*ndcol*L_ENTRY*sizeof(devPtrLx[0]),
                              cudaMemcpyHostToDevice,
                              Common->gpuStream[gpuid][iDevBuff] );

  if ( cudaErr ) {
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"cudaMemcpyAsync H-D");
    return (0);
  }



  /* make the current stream wait for kernels in previous streams */
  cudaStreamWaitEvent ( Common->gpuStream[gpuid][iDevBuff],
                        Common->updateCKernelsComplete[gpuid], 0 ) ;



  /* create relative map for the descendant */
  createRelativeMapOnDevice ( (Int *)(gpu_p->d_Map_root[gpuid]),
                              (Int *)(gpu_p->d_Ls_root[gpuid]),
                              (Int *)(gpu_p->d_RelativeMap_root[gpuid]),
                               pdi1, 
                               ndrow2,
                               &(Common->gpuStream[gpuid][iDevBuff]) );

  cudaErr = cudaGetLastError();
  if (cudaErr) { 
    CHOLMOD_HANDLE_CUDA_ERROR(cudaErr,"createRelativeMapOnDevice");
  }



  /* set cuBlas stream  */
  cublasStatus = cublasSetStream (Common->cublasHandle[gpuid], Common->gpuStream[gpuid][iDevBuff]) ;
  if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
    ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS stream") ;
    return(0);
  }




  /* 
   * Perform DSYRK on GPU for current descendant
   */

  alpha  = 1.0 ;
  beta   = 0.0 ;

#ifdef REAL
  cublasStatus = cublasDsyrk (Common->cublasHandle[gpuid],
                              CUBLAS_FILL_MODE_LOWER,
                              CUBLAS_OP_N,
                              (int) ndrow1,
                              (int) ndcol,    				/* N, K: L1 is ndrow1-by-ndcol */
                              &alpha,         				/* ALPHA:  1 */
                              devPtrLx,
                              ndrow2,         				/* A, LDA: L1, ndrow2 */
                              &beta,          				/* BETA:   0 */
                              devPtrC,
                              ndrow2) ;       				/* C, LDC: C1 */
#else
  cublasStatus = cublasZherk (Common->cublasHandle[gpuid],
        		      CUBLAS_FILL_MODE_LOWER,
  		              CUBLAS_OP_N,
		              (int) ndrow1,
	                      (int) ndcol,    				/* N, K: L1 is ndrow1-by-ndcol*/
		              &alpha,         				/* ALPHA:  1 */
		              (const cuDoubleComplex *) devPtrLx,
		              ndrow2,         				/* A, LDA: L1, ndrow2 */
		              &beta,          				/* BETA:   0 */
		              (cuDoubleComplex *) devPtrC,
		              ndrow2) ;       				/* C, LDC: C1 */
#endif


  if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
    ERROR (CHOLMOD_GPU_PROBLEM, "GPU cublasDsyrk error") ;
    return(0);
  }




  /*
   * Perform DSYRK on GPU for current descendant
   */
  if (ndrow3 > 0)
  {

#ifndef REAL
    cuDoubleComplex calpha  = {1.0,0.0} ;
    cuDoubleComplex cbeta   = {0.0,0.0} ;
#endif

#ifdef REAL
    alpha  = 1.0 ;
    beta   = 0.0 ;

    cublasStatus = cublasDgemm (Common->cublasHandle[gpuid],
                                CUBLAS_OP_N, CUBLAS_OP_T,
                                ndrow3, ndrow1, ndcol,          	/* M, N, K */
                                &alpha,                         	/* ALPHA:  1 */
                                devPtrLx + L_ENTRY*(ndrow1),    	/* A, LDA: L2*/
                                ndrow2,                         	/* ndrow */
                                devPtrLx,                       	/* B, LDB: L1 */
                                ndrow2,                         	/* ndrow */
                                &beta,                          	/* BETA:   0 */
                                devPtrC + L_ENTRY*ndrow1,       	/* C, LDC: C2 */
                                ndrow2) ;
#else
    cublasStatus = cublasZgemm (Common->cublasHandle[gpuid],
    	 		        CUBLAS_OP_N, CUBLAS_OP_C,
		                ndrow3, ndrow1, ndcol,          	/* M, N, K */
		                &calpha,                        	/* ALPHA:  1 */
		                (const cuDoubleComplex*) devPtrLx + ndrow1,
		                ndrow2,                         	/* ndrow */
		                (const cuDoubleComplex *) devPtrLx,
		                ndrow2,                         	/* ndrow */
		                &cbeta,                         	/* BETA:   0 */
		                (cuDoubleComplex *)devPtrC + ndrow1,
		                ndrow2) ;
#endif



    if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU cublasDgemm error") ;
      return(0);
    }

  }




  /*
   * Assemble the update C on the devicet
   */
#ifdef REAL
  addUpdateOnDevice ( gpu_p->d_A_root[gpuid][0], 
		      devPtrC,
                      gpu_p->d_RelativeMap_root[gpuid], 
		      ndrow1, 
		      ndrow2, 
		      nsrow,
                      &(Common->gpuStream[gpuid][iDevBuff]) );
#else
  addComplexUpdateOnDevice ( gpu_p->d_A_root[gpuid][0], 
			     devPtrC,
        		     gpu_p->d_RelativeMap_root[gpuid], 
			     ndrow1, 
			     ndrow2, 
			     nsrow,
        		     &(Common->gpuStream[gpuid][iDevBuff]) );
#endif



  cudaErr = cudaGetLastError();
  if (cudaErr) { 
    ERROR (CHOLMOD_GPU_PROBLEM,"\naddUpdateOnDevice error!\n"); 
    return (0) ; 
  }



  /* record event indicating that kernels for descendant are complete */
  cudaEventRecord ( Common->updateCKernelsComplete[gpuid], Common->gpuStream[gpuid][iDevBuff]);
  cudaEventRecord ( Common->updateCBuffersFree[gpuid][iHostBuff], Common->gpuStream[gpuid][iDevBuff]);



  return (1) ;
}










/*
 *  Function:
 *    gpu_final_assembly_root
 *
 *  Description:
 *    Sum all schur-comlement updates (computed on the GPU) to the supernode
 */
void TEMPLATE2 (CHOLMOD (gpu_final_assembly_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,
   double *Lx,
   int *iHostBuff,
   int *iDevBuff,   
   Int psx,
   Int nscol,
   Int nsrow,
   int supernodeUsedGPU,
   int gpuid
   )
{
  /* local variables */
  Int iidx, i, j, iHostBuff2, iDevBuff2;
  cudaError_t cudaErr ;
  int numThreads;



  numThreads = Common->ompNumThreads;



  /* only if descendant was assembled on GPU */
  if ( supernodeUsedGPU ) {

    /* set host/device buffer coutners */
    *iHostBuff = (Common->ibuffer[gpuid])%CHOLMOD_HOST_SUPERNODE_BUFFERS;
    *iDevBuff = (Common->ibuffer[gpuid])%2;


    /* only if descendant is large enough for GPU */
    if ( nscol * L_ENTRY >= CHOLMOD_POTRF_LIMIT ) {

      /* wait until a buffer is free */
      cudaEventSynchronize ( Common->updateCBuffersFree[gpuid][*iHostBuff] );
      cudaErr = cudaGetLastError();
      if (cudaErr) { 
        ERROR (CHOLMOD_GPU_PROBLEM,"\nsynchronize error!\n"); 
        return; 
      }


      /* copy update assembled on CPU to a pinned buffer */
#pragma omp parallel for num_threads(numThreads)   \
  private(iidx) if (nscol>32)

      for ( j=0; j<nscol; j++ ) {
        for ( i=j; i<nsrow*L_ENTRY; i++ ) {
          iidx = j*nsrow*L_ENTRY+i;
          gpu_p->h_Lx_root[gpuid][*iHostBuff][iidx] = Lx[psx*L_ENTRY+iidx];
        }
      }


      /* H2D transfer of update assembled on CPU */
      cudaMemcpyAsync ( gpu_p->d_A_root[gpuid][1], gpu_p->h_Lx_root[gpuid][*iHostBuff],
                        nscol*nsrow*L_ENTRY*sizeof(double),
                        cudaMemcpyHostToDevice,
                        Common->gpuStream[gpuid][*iDevBuff] );

      cudaErr = cudaGetLastError();
      if (cudaErr) { 
        ERROR (CHOLMOD_GPU_PROBLEM,"\nmemcopy H-D error!\n"); 
        return; 
      }

    } /* end if descendant large enough */


    /* update buffer counters */
    Common->ibuffer[gpuid]++;
    iHostBuff2 = (Common->ibuffer[gpuid])%CHOLMOD_HOST_SUPERNODE_BUFFERS;
    iDevBuff2 = (Common->ibuffer[gpuid])%2;


    /* wait for all kernels to complete */
    cudaEventSynchronize( Common->updateCKernelsComplete[gpuid] );


    /* copy assembled Schur-complement updates computed on GPU */
    cudaMemcpyAsync ( gpu_p->h_Lx_root[gpuid][iHostBuff2], 
		      gpu_p->d_A_root[gpuid][0],
                      nscol*nsrow*L_ENTRY*sizeof(double),
                      cudaMemcpyDeviceToHost,
                      Common->gpuStream[gpuid][iDevBuff2] );

    cudaErr = cudaGetLastError();
    if (cudaErr) { 
      ERROR (CHOLMOD_GPU_PROBLEM,"\nmemcopy D-H error!\n");
      return ; 
    }


    /* need both H2D and D2H copies to be complete */
    cudaDeviceSynchronize();



    /* only if descendant large enough for GPU */
    if ( nscol * L_ENTRY >= CHOLMOD_POTRF_LIMIT ) {

      /* 
       * sum updates from cpu and device on device 
       */
#ifdef REAL
      sumAOnDevice ( gpu_p->d_A_root[gpuid][1], gpu_p->d_A_root[gpuid][0], -1.0, nsrow, nscol );
#else
      sumComplexAOnDevice ( gpu_p->d_A_root[gpuid][1], gpu_p->d_A_root[gpuid][0], -1.0, nsrow, nscol );
#endif


      cudaErr = cudaGetLastError();
      if (cudaErr) { 
        ERROR (CHOLMOD_GPU_PROBLEM,"\nsumAonDevice error!\n");
        return; 
      }


      /* place final assembled supernode in pinned buffer */
      #pragma omp parallel for num_threads(numThreads) private(iidx) if (nscol>32)
      for ( j=0; j<nscol; j++ ) {
        for ( i=j*L_ENTRY; i<nscol*L_ENTRY; i++ ) {
          iidx = j*nsrow*L_ENTRY+i;
          gpu_p->h_Lx_root[gpuid][*iHostBuff][iidx] -=
            gpu_p->h_Lx_root[gpuid][iHostBuff2][iidx];
        }
      }

    } /* end if descendant large enough */
    /* if descendant too small assemble on CPU */
    else
    {

      /* assemble with CPU updates */
      #pragma omp parallel for num_threads(numThreads) private(iidx) if (nscol>32)
      for ( j=0; j<nscol; j++ ) {
        for ( i=j*L_ENTRY; i<nsrow*L_ENTRY; i++ ) {
          iidx = j*nsrow*L_ENTRY+i;
          Lx[psx*L_ENTRY+iidx] -= gpu_p->h_Lx_root[gpuid][iHostBuff2][iidx];
        }
      }

    }


  } /* end if descendant assembled on GPU */


}










/*
 *  Function:
 *    gpu_lower_potrf_root
 *
 *  Description:
 *    computes cholesky factoriation for a supernode
 *    Performs one task:
 *      1. DPOTRF
 */
int TEMPLATE2 (CHOLMOD (gpu_lower_potrf_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,   
   double *Lx,     				
   Int *info,	   				
   Int nscol2,     				
   Int nsrow,      				
   Int psx,        				
   int gpuid
   )
{
  /* local variables */
  int ilda, ijb, iinfo = 0;
  Int j, n, jb, nb, nsrow2, gpu_lda, lda, gpu_ldb;
  double alpha, beta;
  double *devPtrA, *devPtrB, *A;
  cudaError_t cudaErr ;
  cublasStatus_t cublasStatus ;





  /* early exit if descendant is too small for cuBlas */
  if (nscol2 * L_ENTRY < CHOLMOD_POTRF_LIMIT)
  {
    return (0) ;
  }

  
  /* set dimnsions & strides */
  nsrow2 = nsrow - nscol2 ;
  n  = nscol2 ;
  gpu_lda = ((nscol2+31)/32)*32 ;
  lda = nsrow ;
  gpu_ldb = 0 ;
  if (nsrow2 > 0) {
    gpu_ldb = ((nsrow2+31)/32)*32 ;
  }


  /* heuristic to get the block size depending of the problem size */
  nb = 128 ;
  if (nscol2 > 4096) nb = 256 ;
  if (nscol2 > 8192) nb = 384 ;


  /* set device pointers */
  A = gpu_p->h_Lx_root[gpuid][(Common->ibuffer[gpuid]+CHOLMOD_HOST_SUPERNODE_BUFFERS-1)%CHOLMOD_HOST_SUPERNODE_BUFFERS];
  devPtrA = gpu_p->d_Lx_root[gpuid][0];
  devPtrB = gpu_p->d_Lx_root[gpuid][1];



  /* copy A from device to device */
  cudaErr = cudaMemcpy2DAsync ( devPtrA,
                                gpu_lda * L_ENTRY * sizeof (devPtrA[0]),
                                gpu_p->d_A_root[gpuid][1],
                                nsrow * L_ENTRY * sizeof (Lx[0]),
                                nscol2 * L_ENTRY * sizeof (devPtrA[0]),
                                nscol2,
                                cudaMemcpyDeviceToDevice,
                                Common->gpuStream[gpuid][0] );

  if ( cudaErr ){
    ERROR ( CHOLMOD_GPU_PROBLEM, "GPU memcopy device to device");
  }



  /* copy B in advance, for gpu_triangular_solve */
  if (nsrow2 > 0) {
    cudaErr = cudaMemcpy2DAsync (devPtrB,
                                 gpu_ldb * L_ENTRY * sizeof (devPtrB [0]),
                                 gpu_p->d_A_root[gpuid][1] + L_ENTRY*nscol2,
                                 nsrow * L_ENTRY * sizeof (Lx [0]),
                                 nsrow2 * L_ENTRY * sizeof (devPtrB [0]),
                                 nscol2,
                                 cudaMemcpyDeviceToDevice,
                                 Common->gpuStream[gpuid][0]) ;
    if (cudaErr) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
    }
  }



  /* define the dpotrf stream */
  cublasStatus = cublasSetStream (Common->cublasHandle[gpuid], Common->gpuStream[gpuid][0]) ;
  if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
    ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS stream") ;
  }



  /* 
   * block Cholesky factorization of S 
   */
  /* loop over blocks */
  for (j = 0 ; j < n ; j += nb) {

    jb = nb < (n-j) ? nb : (n-j) ;     


    /* 
     * Perform DSYRK on GPU 
     */
    alpha = -1.0;
    beta  = 1.0;

#ifdef REAL
    cublasStatus = cublasDsyrk (Common->cublasHandle[gpuid],
                                CUBLAS_FILL_MODE_LOWER, 
 		                CUBLAS_OP_N, 
			        jb, 
			        j,
                                &alpha, devPtrA + j, 
				gpu_lda,
                                &beta,  devPtrA + j + j*gpu_lda, 
				gpu_lda);
#else
    cublasStatus = cublasZherk (Common->cublasHandle[gpuid],
    			        CUBLAS_FILL_MODE_LOWER, 
			 	CUBLAS_OP_N, 
				jb, 
				j,
   			        &alpha, (cuDoubleComplex*)devPtrA + j,
		                gpu_lda,
		                &beta,
		                (cuDoubleComplex*)devPtrA + j + j*gpu_lda,
		                gpu_lda) ;
#endif

    if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU cublasDsyrk error") ;
    }


    /* record end of cublasDsyrk */
    cudaErr = cudaEventRecord (Common->cublasEventPotrf[gpuid][0], Common->gpuStream[gpuid][0]) ;
    if (cudaErr) {
      ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event failure") ;
    }


    /* wait for cublasDsyrk to end */
    cudaErr = cudaStreamWaitEvent (Common->gpuStream[gpuid][1], Common->cublasEventPotrf[gpuid][0], 0) ;
    if (cudaErr) {
      ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event failure") ;
    }


    /* copy back the jb columns on two different streams */
    cudaErr = cudaMemcpy2DAsync (A + L_ENTRY*(j + j*lda),
                                 lda * L_ENTRY * sizeof (double),
                                 devPtrA + L_ENTRY*(j + j*gpu_lda),
                                 gpu_lda * L_ENTRY * sizeof (double),
                                 L_ENTRY * sizeof (double)*jb,
                                 jb,
                                 cudaMemcpyDeviceToHost,
                                 Common->gpuStream[gpuid][1]) ;

    if (cudaErr) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy from device") ;
    }



    /*
     * Perform DGEMM on GPU   
     */
    if ((j+jb) < n) {


#ifdef REAL
      alpha = -1.0 ;
      beta  = 1.0 ;
      cublasStatus = cublasDgemm (Common->cublasHandle[gpuid],
                                  CUBLAS_OP_N, 
				  CUBLAS_OP_T,
                                  (n-j-jb), 
				  jb, 
				  j,
                                  &alpha,
                                  devPtrA + (j+jb), 
				  gpu_lda,
                                  devPtrA + (j), 
				  gpu_lda,
                                  &beta,
                                  devPtrA + (j+jb + j*gpu_lda), 
 		                  gpu_lda);
#else
      cuDoubleComplex calpha = {-1.0,0.0} ;
      cuDoubleComplex cbeta  = { 1.0,0.0} ;

      cublasStatus = cublasZgemm (Common->cublasHandle[gpuid],
	                   	  CUBLAS_OP_N, 
				  CUBLAS_OP_C,
   		                  (n-j-jb), 
				  jb, 
				  j,
		                  &calpha,
		                  (cuDoubleComplex*)devPtrA + (j+jb),
		                  gpu_lda,
		                  (cuDoubleComplex*)devPtrA + (j),
		                  gpu_lda,
		                  &cbeta,
		                  (cuDoubleComplex*)devPtrA + (j+jb + j*gpu_lda),
		                  gpu_lda ) ;
#endif


      if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
      }

    }


    /* synchronize stream */
    cudaErr = cudaStreamSynchronize (Common->gpuStream[gpuid][1]) ;
    if (cudaErr) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
    }


    /* compute the Cholesky factorization of the jbxjb block on the CPU */
    ilda = (int) lda ;
    ijb  = jb ;

#ifdef REAL
    LAPACK_DPOTRF ("L", &ijb, A + L_ENTRY * (j + j*lda), &ilda, &iinfo) ;
#else
    LAPACK_ZPOTRF ("L", &ijb, A + L_ENTRY * (j + j*lda), &ilda, &iinfo) ;
#endif


    /* get parameter that determines if it is positive definite */
    *info = iinfo ;
    if (*info != 0) {
      *info = *info + j ;
      break ;
    }



    /* copy the result back to the GPU */
    cudaErr = cudaMemcpy2DAsync (devPtrA + L_ENTRY*(j + j*gpu_lda),
                                 gpu_lda * L_ENTRY * sizeof (double),
                                 A + L_ENTRY * (j + j*lda),
                                 lda * L_ENTRY * sizeof (double),
                                 L_ENTRY * sizeof (double) * jb,
                                 jb,
                                 cudaMemcpyHostToDevice,
                                 Common->gpuStream[gpuid][0]) ;

    if (cudaErr) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
    }



    /* 
     * Perform DTRSM on GPU
     */ 
    if ((j+jb) < n) {

#ifdef REAL
      alpha  = 1.0 ;
      cublasStatus = cublasDtrsm (Common->cublasHandle[gpuid],
                                  CUBLAS_SIDE_RIGHT,
                                  CUBLAS_FILL_MODE_LOWER,
                                  CUBLAS_OP_T, 
				  CUBLAS_DIAG_NON_UNIT,
                                  (n-j-jb), 
				  jb,
                                  &alpha,
                                  devPtrA + (j + j*gpu_lda), 
				  gpu_lda,
                                  devPtrA + (j+jb + j*gpu_lda), 
				  gpu_lda);
#else
      cuDoubleComplex calpha  = {1.0,0.0};

      cublasStatus = cublasZtrsm (Common->cublasHandle[gpuid],
		                  CUBLAS_SIDE_RIGHT,
               			  CUBLAS_FILL_MODE_LOWER,
		                  CUBLAS_OP_C, 
				  CUBLAS_DIAG_NON_UNIT,
		                  (n-j-jb), 
				  jb,
		                  &calpha,
		                  (cuDoubleComplex *)devPtrA + (j + j*gpu_lda),
		                  gpu_lda,
		                  (cuDoubleComplex *)devPtrA + (j+jb + j*gpu_lda),
		                  gpu_lda) ;
#endif


      if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
      }


      /* record end of DTRSM  */
      cudaErr = cudaEventRecord (Common->cublasEventPotrf[gpuid][2], Common->gpuStream[gpuid][0]) ;
      if (cudaErr) {
        ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event failure") ;
      }


      /* wait fo end of DTRSM */
      cudaErr = cudaStreamWaitEvent (Common->gpuStream[gpuid][1], Common->cublasEventPotrf[gpuid][2], 0) ;
      if (cudaErr) {
        ERROR (CHOLMOD_GPU_PROBLEM, "CUDA event failure") ;
      }


      /* Copy factored column back to host. */
      cudaErr = cudaMemcpy2DAsync (A + L_ENTRY*(j + jb + j * lda),
                                   lda * L_ENTRY * sizeof (double),
                                   devPtrA + L_ENTRY*
                                   (j + jb + j * gpu_lda),
                                   gpu_lda * L_ENTRY * sizeof (double),
                                   L_ENTRY * sizeof (double)*
                                   (n - j - jb), jb,
                                   cudaMemcpyDeviceToHost,
                                   Common->gpuStream[gpuid][1]) ;

      if (cudaErr) {
        ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy to device") ;
      }
    }

  } /* end loop over blocks */



  return (1) ;
}










/*
 *  Function:
 *    gpu_triangular_solve_root
 *
 *  Description:
 *    Computes triangular solve for a supernode.
 *    Performs one task:
 *      1. DTRSM
 *
 */
int TEMPLATE2 (CHOLMOD (gpu_triangular_solve_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,
   double *Lx,
   Int nsrow2,     
   Int nscol2,     
   Int nsrow,      
   Int psx,        
   int gpuid
   )
{
  /* local variables */ 
  int i, j, iwrap, ibuf = 0, iblock = 0, iHostBuff, numThreads;
  Int iidx, gpu_lda, gpu_ldb, gpu_rowstep, gpu_row_max_chunk, gpu_row_chunk, gpu_row_start = 0;
  double *devPtrA, *devPtrB;
  cudaError_t cudaErr;
  cublasStatus_t cublasStatus;




  numThreads = Common->ompNumThreads;

 

  /* early exit */
  if ( nsrow2 <= 0 )
  {
    return (0) ;
  }


  /* initialize parameters */
  iHostBuff = (Common->ibuffer[gpuid]+CHOLMOD_HOST_SUPERNODE_BUFFERS-1) % CHOLMOD_HOST_SUPERNODE_BUFFERS;
  gpu_lda = ((nscol2+31)/32)*32 ;
  gpu_ldb = ((nsrow2+31)/32)*32 ;
#ifdef REAL
    double alpha  = 1.0 ;
    gpu_row_max_chunk = 768;
#else
    cuDoubleComplex calpha  = {1.0,0.0} ;
    gpu_row_max_chunk = 256;
#endif


  /* initialize device pointers */
  devPtrA = gpu_p->d_Lx_root[gpuid][0];
  devPtrB = gpu_p->d_Lx_root[gpuid][1];


  /* make sure the copy of B has completed */
  cudaStreamSynchronize( Common->gpuStream[gpuid][0] );




  /*
   * Perform blocked TRSM:
   * 1. compute DTRSM
   * 2. copy Lx from pinned to host memory
   * Hide copies behind compute.
   */
  /* loop over blocks */
  while ( gpu_row_start < nsrow2 ) {

    /* set block dimensions */
    gpu_row_chunk = nsrow2 - gpu_row_start;
    if ( gpu_row_chunk  > gpu_row_max_chunk ) {
      gpu_row_chunk = gpu_row_max_chunk;
    }


    cublasStatus = cublasSetStream ( Common->cublasHandle[gpuid], Common->gpuStream[gpuid][ibuf] );

    if ( cublasStatus != CUBLAS_STATUS_SUCCESS ) {
      ERROR ( CHOLMOD_GPU_PROBLEM, "GPU CUBLAS stream");
    }


    /*
     * Perform DTRSM on GPU
     */
#ifdef REAL
    cublasStatus = cublasDtrsm (Common->cublasHandle[gpuid],
                                CUBLAS_SIDE_RIGHT,
                                CUBLAS_FILL_MODE_LOWER,
                                CUBLAS_OP_T,
                                CUBLAS_DIAG_NON_UNIT,
                                gpu_row_chunk,
                                nscol2,
                                &alpha,
                                devPtrA,
                                gpu_lda,
                                devPtrB + gpu_row_start,
                                gpu_ldb) ;
#else
    cublasStatus = cublasZtrsm (Common->cublasHandle[gpuid],
                                CUBLAS_SIDE_RIGHT,
                                CUBLAS_FILL_MODE_LOWER,
                                CUBLAS_OP_C,
                                CUBLAS_DIAG_NON_UNIT,
                                gpu_row_chunk,
                                nscol2,
                                &calpha,
                                (const cuDoubleComplex *) devPtrA,
                                gpu_lda,
                                (cuDoubleComplex *)devPtrB + gpu_row_start ,
                                gpu_ldb) ;
#endif


    if (cublasStatus != CUBLAS_STATUS_SUCCESS) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU CUBLAS routine failure") ;
    }


    /* copy result back to the CPU */
    cudaErr = cudaMemcpy2DAsync ( gpu_p->h_Lx_root[gpuid][iHostBuff] +
                                  L_ENTRY*(nscol2+gpu_row_start),
                                  nsrow * L_ENTRY * sizeof (Lx [0]),
                                  devPtrB + L_ENTRY*gpu_row_start,
                                  gpu_ldb * L_ENTRY * sizeof (devPtrB [0]),
                                  gpu_row_chunk * L_ENTRY *
                                  sizeof (devPtrB [0]),
                                  nscol2,
                                  cudaMemcpyDeviceToHost,
                                  Common->gpuStream[gpuid][ibuf]);

    if (cudaErr) {
      ERROR (CHOLMOD_GPU_PROBLEM, "GPU memcopy from device") ;
    }


    /* record end of copy */
    cudaEventRecord ( Common->updateCBuffersFree[gpuid][ibuf], Common->gpuStream[gpuid][ibuf] );


    /* update block dimensions */
    gpu_row_start += gpu_row_chunk;
    ibuf++;
    ibuf = ibuf % CHOLMOD_HOST_SUPERNODE_BUFFERS;
    iblock ++;


    /* only if enough available host buffers */
    if ( iblock >= CHOLMOD_HOST_SUPERNODE_BUFFERS ) {

      Int gpu_row_start2 ;
      Int gpu_row_end ;


      cudaErr = cudaEventSynchronize ( Common->updateCBuffersFree[gpuid][iblock%CHOLMOD_HOST_SUPERNODE_BUFFERS] );
      if ( cudaErr ) {
        printf ("ERROR cudaEventSynchronize\n");
      }


      gpu_row_start2 = nscol2 + (iblock-CHOLMOD_HOST_SUPERNODE_BUFFERS)*gpu_row_max_chunk;
      gpu_row_end = gpu_row_start2+gpu_row_max_chunk;

      if ( gpu_row_end > nsrow ) gpu_row_end = nsrow;


      /* copy Lx from pinned to host memory */
      #pragma omp parallel for num_threads(numThreads) private(iidx) if ( nscol2 > 32 )
      for ( j=0; j<nscol2; j++ ) {
        for ( i=gpu_row_start2*L_ENTRY; i<gpu_row_end*L_ENTRY; i++ ) {
          iidx = j*nsrow*L_ENTRY+i;
          Lx[psx*L_ENTRY+iidx] = gpu_p->h_Lx_root[gpuid][iHostBuff][iidx];
        }
      }

    } /* end if enough buffers */

  } /* end while loop */




  /* Convenient to copy the L1 block here */
  #pragma omp parallel for num_threads(numThreads) private ( iidx ) if ( nscol2 > 32 )
  for ( j=0; j<nscol2; j++ ) {
    for ( i=j*L_ENTRY; i<nscol2*L_ENTRY; i++ ) {
      iidx = j*nsrow*L_ENTRY + i;
      Lx[psx*L_ENTRY+iidx] = gpu_p->h_Lx_root[gpuid][iHostBuff][iidx];
    }
  }




  /* now account for the last HSTREAMS buffers */
  for ( iwrap=0; iwrap<CHOLMOD_HOST_SUPERNODE_BUFFERS; iwrap++ ) {

    int i, j;
    Int gpu_row_start2 = nscol2 + (iblock-CHOLMOD_HOST_SUPERNODE_BUFFERS)*gpu_row_max_chunk;

    if (iblock-CHOLMOD_HOST_SUPERNODE_BUFFERS >= 0 && gpu_row_start2 < nsrow ) {

      Int iidx;
      Int gpu_row_end = gpu_row_start2+gpu_row_max_chunk;
      if ( gpu_row_end > nsrow ) gpu_row_end = nsrow;

      cudaEventSynchronize ( Common->updateCBuffersFree[gpuid][iblock%CHOLMOD_HOST_SUPERNODE_BUFFERS] );
         
 
      /* copy Lx from pinned to host memory */
      #pragma omp parallel for num_threads(numThreads) private(iidx) if ( nscol2 > 32 )
      for ( j=0; j<nscol2; j++ ) {
        for ( i=gpu_row_start2*L_ENTRY; i<gpu_row_end*L_ENTRY; i++ ) {
          iidx = j*nsrow*L_ENTRY+i;
          Lx[psx*L_ENTRY+iidx] = gpu_p->h_Lx_root[gpuid][iHostBuff][iidx];
        }
      }
    }
  iblock++;
  } /* end loop over wrap */



  return (1) ;
}











/*
 *  Function:
 *    gpu_copy_supernode_root
 *
 *  Description:
 *    Copies supernode from pinned to host memory.
 *    Case triangular_solve is not called..
 */
void TEMPLATE2 (CHOLMOD (gpu_copy_supernode_root))
  (
   cholmod_common *Common,
   cholmod_gpu_pointers *gpu_p,
   double *Lx,
   Int psx,
   Int nscol,
   Int nscol2,
   Int nsrow,
   int supernodeUsedGPU,
   int iHostBuff,
   int gpuid
   )
{
  /* local variables */
  Int iidx, i, j;
  int numThreads;


  numThreads = Common->ompNumThreads;


  /* if supernode large enough for GPU */
  if ( supernodeUsedGPU && nscol2 * L_ENTRY >= CHOLMOD_POTRF_LIMIT ) {
  
    /* synchronize device */
    cudaDeviceSynchronize();

    /* copy Lx from pinned to host memory */
    #pragma omp parallel for num_threads(numThreads) private(iidx,i,j) if (nscol>32)
    for ( j=0; j<nscol; j++ ) {
      for ( i=j*L_ENTRY; i<nscol*L_ENTRY; i++ ) {
        iidx = j*nsrow*L_ENTRY+i;
        Lx[psx*L_ENTRY+iidx] = gpu_p->h_Lx_root[gpuid][iHostBuff][iidx];
      }
    }

  } /* end if statement */


}











#undef REAL
#undef COMPLEX
#undef ZOMPLEX












