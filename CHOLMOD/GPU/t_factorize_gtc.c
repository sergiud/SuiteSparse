/* ========================================================================== */
/* === GPU/t_factorize_gtc.c ================================================ */
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
 *   t_factorize_gtc
 *
 * Description:
 *   Contains functions for factorization
 *   of the GTC algorithm.
 *
 */

/*

  PARALLEL TODO

  1. define a Map for each GPU - need to know the max size reqd. for Map.

  2. remove ndescendants argument.  I don't see that it is needed.

  3. 

 */


/* includes */
#include <string.h>
#include <time.h>
#ifdef MKLROOT
#include <mkl.h>
#endif

/* nvtx markers (for profiling) */
/*#define USE_NVTX*/
#ifdef USE_NVTX
#include "nvToolsExt.h"
#endif

/* undef macros */
#undef L_ENTRY
#undef L_CLEAR
#undef L_ASSIGN
#undef L_MULTADD
#undef L_ASSEMBLE
#undef L_ASSEMBLESUB


/* macros */
#ifdef REAL

/* A, F, and L are all real */
#define L_ENTRY 1
#define L_CLEAR(Lx,p)                   Lx [p] = 0
#define L_ASSIGN(Lx,q, Ax,Az,p)         Lx [q] = Ax [p]
#define L_MULTADD(Lx,q, Ax,Az,p, f)     Lx [q] += Ax [p] * f [0]
#define L_ASSEMBLE(Lx,q,b)              Lx [q] += b [0]
#define L_ASSEMBLESUB(Lx,q,C,p)         Lx [q] -= C [p]

#else

/* A and F are complex or zomplex, L and C are complex */
#define L_ENTRY 2
#define L_CLEAR(Lx,p)                   Lx [2*(p)] = 0 ; Lx [2*(p)+1] = 0
#define L_ASSEMBLE(Lx,q,b)              Lx [2*(q)] += b [0] ;
#define L_ASSEMBLESUB(Lx,q,C,p)         Lx [2*(q)  ] -= C [2*(p)  ] ;           \
                                        Lx [2*(q)+1] -= C [2*(p)+1] ;

#ifdef COMPLEX

/* A, F, L, and C are all complex */
#define L_ASSIGN(Lx,q, Ax,Az,p)         Lx [2*(q)  ] = Ax [2*(p)  ] ;           \
                                        Lx [2*(q)+1] = Ax [2*(p)+1]

#define L_MULTADD(Lx,q, Ax,Az,p, f)     Lx [2*(q)  ] += Ax [2*(p)  ] * f [0] - Ax [2*(p)+1] * f [1] ;           \
                                        Lx [2*(q)+1] += Ax [2*(p)+1] * f [0] + Ax [2*(p)  ] * f [1]

#else

/* A and F are zomplex, L and C is complex */
#define L_ASSIGN(Lx,q, Ax,Az,p)         Lx [2*(q)  ] = Ax [p] ;                 \
                                        Lx [2*(q)+1] = Az [p] ;

#define L_MULTADD(Lx,q, Ax,Az,p, f)     Lx [2*(q)  ] += Ax [p] * f [0] - Az [p] * f [1] ;   \
                                        Lx [2*(q)+1] += Az [p] * f [0] + Ax [p] * f [1]

#endif
#endif











/*
 * Function:
 *   gpu_factorize_gtc
 *
 * Description:
 *   Factorizes top-of-tree branch of elimination tree, where
 *   the branch does not fit the GPU. Utilizes a hybrid algorithm
 *   presented at GTC14.
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_factorize_gtc))
  (
    cholmod_common *Common,
    cholmod_factor *L,
    cholmod_gpu_pointers *gpu_p,
    Int branch,
    Int nsuper,
    Int n,
    Int stype,
    Int Apacked,
    Int Fpacked,
    Int *Lpos,
    Int *Lpos_save,
    Int *Lpi,
    Int *Lpx,
    Int *Ls,
    Int *Ap,
    Int *Ai,
    Int *Anz,
    Int *Fp,
    Int *Fi,
    Int *Fnz,
    Int *Super,
    Int *SuperMap,
    Int *Map,
    Int *RelativeMap,
    Int *ndescendants,
    Int *Previous,
    Int *Next,
    Int *Next_save,
    Int *Head,
    Int *supernode_levels,
    Int *supernode_levels_ptrs,
    Int *supernode_levels_branch_ptrs,
    Int *supernode_num_levels,
    double *beta,
    double *C,
    double *Lx,
    double *Ax,
    double *Az,
    double *Fx,
    double *Fz
   )
{
  
  
  Int repeat_supernode = FALSE ;

  Int gpuid, nscol3;
 
  /* set device 0 */
  gpu_p->gpuid = 0;
  gpuid = gpu_p->gpuid;
  nscol3 = 0;
  /* initialize GPU */
  TEMPLATE2 (CHOLMOD (gpu_init_gtc))(Common, gpu_p, L, Lpi, nsuper, n,gpuid);

  /* set device 1 */
  gpu_p->gpuid = 1;
  gpuid = gpu_p->gpuid;
  nscol3 = 0;
  /* initialize GPU */
  TEMPLATE2 (CHOLMOD (gpu_init_gtc))(Common, gpu_p, L, Lpi, nsuper, n,gpuid);

  /* set device 2 */
  gpu_p->gpuid = 2;
  gpuid = gpu_p->gpuid;
  nscol3 = 0;
  /* initialize GPU */
  TEMPLATE2 (CHOLMOD (gpu_init_gtc))(Common, gpu_p, L, Lpi, nsuper, n,gpuid);

  /* set device 3 */
  gpu_p->gpuid = 3;
  gpuid = gpu_p->gpuid;
  nscol3 = 0;
  /* initialize GPU */
  TEMPLATE2 (CHOLMOD (gpu_init_gtc))(Common, gpu_p, L, Lpi, nsuper, n,gpuid);

  /* set device 0 */
  gpu_p->gpuid = 1;
  gpuid = gpu_p->gpuid;

  const double one[2] = {1.0, 0.0};
  const double zero[2] = {0.0, 0.0};



#define THREADS_OUTER 4
#define THREADS_INNER 1  
  omp_set_nested(TRUE);


  Int level, start, end, node;

  /* loop over levels in branch */
  for(level = 0; level < supernode_num_levels[branch]; level++)
  {  

    start = supernode_levels_ptrs[supernode_levels_branch_ptrs[branch]+level];
    end = supernode_levels_ptrs[supernode_levels_branch_ptrs[branch]+level+1];


    /* loop over supernodes */
#pragma omp parallel for num_threads(THREADS_OUTER)
    for(node = start; node < end; node++)
    {

      /* local variables */
	int i, j, k, gpuid;
      Int px, pk, pf, p, q, d, s, ss, ndrow, ndrow1, ndrow2, ndrow3, ndcol, nsrow, nsrow2, nscol, nscol2, nscol3, 
	kd1, kd2, k1, k2, psx, psi, pdx, pdx1, pdi, pdi1, pdi2, pdend, psend, pfend, pend, dancestor, sparent, imap, 
	idescendant, dnext, dlarge, iHostBuff, iDevBuff, skips, skip_max,dsmall, tail, info,
	GPUavailable, mapCreatedOnGpu, supernodeUsedGPU;
      Int ndescendants;
      nscol3 = 0;
      cudaError_t cuErr;
      
      Int thd_d[CHOLMOD_OMP_NUM_THREADS];
      Int thd_kd1[CHOLMOD_OMP_NUM_THREADS];
      Int thd_kd2[CHOLMOD_OMP_NUM_THREADS];
      Int thd_ndcol[CHOLMOD_OMP_NUM_THREADS];
      Int thd_pdi[CHOLMOD_OMP_NUM_THREADS];
      Int thd_pdx[CHOLMOD_OMP_NUM_THREADS];
      Int thd_pdend[CHOLMOD_OMP_NUM_THREADS];
      Int thd_ndrow[CHOLMOD_OMP_NUM_THREADS];
      Int thd_p[CHOLMOD_OMP_NUM_THREADS];
      Int thd_pdi1[CHOLMOD_OMP_NUM_THREADS];
      Int thd_pdi2[CHOLMOD_OMP_NUM_THREADS];
      Int thd_pdx1[CHOLMOD_OMP_NUM_THREADS];
      Int thd_ndrow1[CHOLMOD_OMP_NUM_THREADS];
      Int thd_ndrow2[CHOLMOD_OMP_NUM_THREADS];
      Int thd_ndrow3[CHOLMOD_OMP_NUM_THREADS];
      double *thd_C[CHOLMOD_OMP_NUM_THREADS];
  
      


      /* obtain thread id */
      int outer_tid = 0;
      outer_tid = omp_get_thread_num();

      gpuid = outer_tid;
      
      cudaSetDevice(gpuid);


      /* get supernode dimensions */
      s = supernode_levels[node];
      k1 = Super [s];            		/* s contains columns k1 to k2-1 of L */
      k2 = Super [s+1];
      nscol = k2 - k1;           		/* # of columns in all of s */
      psi = Lpi [s];             		/* pointer to first row of s in Ls */
      psx = Lpx [s];             		/* pointer to first row of s in Lx */
      psend = Lpi [s+1];         		/* pointer just past last row of s in Ls */
      nsrow = psend - psi;       		/* # of rows in all of s */
      pend = psx + nsrow * nscol;       	/* s is nsrow-by-nscol */
      pk = psx;

      Int Map[20000000];
      Int RelativeMap[20000000];
      double *rootC = (double*) malloc (L->maxcsize*sizeof(double));

      /* construct the scattered Map for supernode s */
#pragma omp parallel for num_threads(THREADS_INNER) if ( nsrow > 128 )
      for (k = 0 ; k < nsrow ; k++)
        Map [Ls [psi + k]] = k ;
  



      /* reorder descendants in supernode by descreasing size */
      TEMPLATE2 (CHOLMOD (gpu_reorder_descendants_gtc))(Common, gpu_p, Lpi, Lpos, Super, Head, &tail, Next, Previous,
                                                        &ndescendants/*[0]*/, &mapCreatedOnGpu, s, gpuid );




      /* copy matrix into supernode s (lower triangular part only) */
#pragma omp parallel for private ( p, pend, pfend, pf, i, j, imap, q ) num_threads(THREADS_INNER) if ( k2-k1 > 64 )
      for (k = k1 ; k < k2 ; k++) 
      {
        /* copy the kth column of A into the supernode */
        if (stype != 0) 
        {
          p = Ap [k] ;
          pend = (Apacked) ? (Ap [k+1]) : (p + Anz [k]) ;

          for ( ; p < pend ; p++) 
          {
            i = Ai [p] ;
            if (i >= k) 
            {
              imap = Map [i] ;							/* row i of L is located in row Map [i] of s */
              if (imap >= 0 && imap < nsrow) 
              {
                L_ASSIGN (Lx,(imap+(psx+(k-k1)*nsrow)), Ax,Az,p) ;		/* Lx [Map [i] + pk] = Ax [p] ; */
              }
            }
          }
        }
        /* copy the kth column of A*F into the supernode */
        else				
        {
          double fjk[2];
          pf = Fp [k] ;
          pfend = (Fpacked) ? (Fp [k+1]) : (p + Fnz [k]) ;
          for ( ; pf < pfend ; pf++) 
          {
            j = Fi [pf] ;
            L_ASSIGN (fjk,0, Fx,Fz,pf) ;					/* fjk = Fx [pf] ; */
            p = Ap [j] ;
            pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
            for ( ; p < pend ; p++) 
            {
              i = Ai [p] ;
              if (i >= k)
              {
                imap = Map [i] ;
                if (imap >= 0 && imap < nsrow)
                {
                  L_MULTADD (Lx,(imap+(psx+(k-k1)*nsrow)),Ax,Az,p, fjk) ;	/* Lx [Map [i] + pk] += Ax [p] * fjk ; */
                }
              }
            }
          }
        }
      }
  


      /* add beta (only real part) to the diagonal of the supernode, if nonzero */
      if (beta [0] != 0.0)
      {
        pk = psx ;
        for (k = k1 ; k < k2 ; k++)
        {
          L_ASSEMBLE (Lx,pk, beta) ;	/* Lx [pk] += beta [0] ; */
          pk += nsrow + 1 ;       	/* advance to the next diagonal entry */
        }
      }







      /* save/restore the list of supernodes */
      if (!repeat_supernode)
      {
        for (d = Head [s] ; d != EMPTY ; d = Next [d])
        {
          Lpos_save [d] = Lpos [d] ;
          Next_save [d] = Next [d] ;
        }
      }
      else
      {
        for (d = Head [s] ; d != EMPTY ; d = Next [d])
        {
          Lpos [d] = Lpos_save [d] ;
          Next [d] = Next_save [d] ;
        }
      }



      /* initialize the buffer counter */
      Common->ibuffer[gpuid] = 0;
      supernodeUsedGPU = 0;
      idescendant = 0;
      d = Head[s];
      dnext = d;
      dlarge = Next[d];
      dsmall = tail;
      GPUavailable = 1;
      skip_max = 0;                               /* initially, while performing CHOLMOD_OMP_NUM_THREADS BLAS 
						     call simultaneously, don't skip any calls to 
						     cudaEventQuery */
      skips = 0;

      int THREADS = THREADS_INNER;  

#ifdef MKLROOT
      mkl_set_num_threads_local(1);     
#else
      openblas_set_num_threads(1);     
#endif




      /* loop over descendants d of supernode s */
      while( (idescendant < ndescendants/*[0]*/) ) {
  
	iHostBuff = (Common->ibuffer[gpuid]) % CHOLMOD_HOST_SUPERNODE_BUFFERS;

        if ( (THREADS > 1) && (ndescendants/*[0]*/-idescendant < CHOLMOD_OMP_NUM_THREADS*4) ) {
          THREADS=1;

#ifdef MKLROOT
          mkl_set_num_threads_local(CHOLMOD_OMP_NUM_THREADS);
#else
          openblas_set_num_threads(CHOLMOD_OMP_NUM_THREADS);
#endif

	  skip_max = CHOLMOD_GPU_SKIP;
        }

        /* get next descendant */  
        if ( idescendant > 0 )  {
          if ( GPUavailable == -1 || skips > 0 ) {       /* -1 indicates that all GPU-eligible large supernodes have been processed 
								so don't need to call expensive cudaEventQuery */
            d = dsmall;
            dsmall = Previous[dsmall];
            (skips)--;
          }
          else {
            cuErr = cudaEventQuery( Common->updateCBuffersFree[gpuid][iHostBuff] );
            if ( cuErr == cudaSuccess ) {
              d = dlarge;
              dlarge = Next[dlarge];
              GPUavailable = 1;                         /* signals that descendant will go to GPU  */
              skips = 0;
            }
            else {                                    
              d = dsmall;
              dsmall = Previous[dsmall];
              GPUavailable = 0;                         /* signals that this descendant will go to the CPU*/
              skips = skip_max;                 /* to reduce the cost of cudaEventQuery - only query after
	                                                  /*processing CHOLMOD_GPU_SKIP descendant supernodes       */
            } 
          }
        }
  
        (idescendant)++;


    
        /* get the size of supernode d */
        kd1 = Super [d] ;      		 	/* d contains cols kd1 to kd2-1 of L */
        kd2 = Super [d+1] ;
        ndcol = kd2 - kd1 ; 		 	/* # of columns in all of d */
        pdi = Lpi [d] ;         		/* pointer to first row of d in Ls */
        pdx = Lpx [d] ;         	 	/* pointer to first row of d in Lx */
        pdend = Lpi [d+1] ;     	 	/* pointer just past last row of d in Ls */
        ndrow = pdend - pdi ;   	 	/* # rows in all of d */
  
        /* find the range of rows of d that affect rows k1 to k2-1 of s */
        p = Lpos [d] ;          	 	/* offset of 1st row of d affecting s */
        pdi1 = pdi + p ;        	 	/* ptr to 1st row of d affecting s in Ls */
        pdx1 = pdx + p ;        	 	/* ptr to 1st row of d affecting s in Lx */
  
        for (pdi2 = pdi1 ; pdi2 < pdend && Ls [pdi2] < k2 ; (pdi2)++) ;
        ndrow1 = pdi2 - pdi1 ;      	 	/* # rows in first part of d */
        ndrow2 = pdend - pdi1 ;     	 	/* # rows in remaining d */
  
        /* construct the update matrix C for this supernode d */
        ndrow3 = ndrow2 - ndrow1 ;  	 	/* number of rows of C2 */
  


      /*
       *  Initialize Supernode
       *
       *  Initializes the supernode with the following steps:
       *
       *  1. clear supernode (Lx) on device
       *  2. create Map for supernode
       *
       */
        if ( GPUavailable == 1) {
          if ( ndrow2 * L_ENTRY >= CHOLMOD_ND_ROW_LIMIT &&
            ndcol * L_ENTRY >= CHOLMOD_ND_COL_LIMIT ) {
              if ( ! mapCreatedOnGpu ) {
                /* initialize supernode (create Map) */
                TEMPLATE2 ( CHOLMOD (gpu_initialize_supernode_gtc))( Common, gpu_p, nscol, nsrow, psi, gpuid );
                mapCreatedOnGpu = 1;
              }
            }
          else {
	    /* intercept the end of GPU-eligible descendant supernodes here */
            GPUavailable = -1;
          }
        }



      /*
       *  Supernode Assembly
       *
       *  Assemble the supernode with the following steps:
       *
       *  1. perform dsyrk
       *  2. perform dgemm
       *  3. perform addUpdate
       *
       */

        if ( GPUavailable!=1 || !TEMPLATE2 (CHOLMOD (gpu_updateC_gtc)) (Common, gpu_p, Lx, C, ndrow1, ndrow2, ndrow, ndcol, nsrow, pdx1, pdi1, gpuid) ) {
	  

	  /* now, actually process up to 32 simultaneous descendent supernodes to better utilize CPU BLAS/cores */
	  
	  int nlocalt = 0;
	  
	  /* make sure we can fit all of our output data in C (likely since C is sized for largest supernode and here we are processing smallest) */
	  double *ptrC = rootC;
	  Int Cused=0;

	  int tid;
	  for ( tid=0; tid<THREADS; tid++ ) {

	    if ( ( idescendant <= ndescendants/*[0]*/ && tid==0) ||
		 ( idescendant < ndescendants/*[0]*/ ) ) {	
	      if ( tid == 0 ) {
		thd_d[tid] = d;
		nlocalt++;
	      }
	      else {
		thd_d[tid] = dsmall;
	      }

	      /* get the size of supernode d */
	      thd_kd1[tid] = Super [thd_d[tid]] ;      		                /* d contains cols kd1 to kd2-1 of L */
	      thd_kd2[tid] = Super [thd_d[tid]+1] ;
	      thd_ndcol[tid] = thd_kd2[tid] - thd_kd1[tid] ; 	                /* # of columns in all of d */
	      thd_pdi[tid] = Lpi [thd_d[tid]] ;         	  	        /* pointer to first row of d in Ls */
	      thd_pdx[tid] = Lpx [thd_d[tid]] ;         	 	        /* pointer to first row of d in Lx */
	      thd_pdend[tid] = Lpi [thd_d[tid]+1] ;     	 	        /* pointer just past last row of d in Ls */
	      thd_ndrow[tid] = thd_pdend[tid] - thd_pdi[tid] ;   	 	/* # rows in all of d */
	      
	      /* find the range of rows of d that affect rows k1 to k2-1 of s */
	      thd_p[tid] = Lpos [thd_d[tid]] ;          	 	        /* offset of 1st row of d affecting s */
	      thd_pdi1[tid] = thd_pdi[tid] + thd_p[tid] ;        	 	/* ptr to 1st row of d affecting s in Ls */
	      thd_pdx1[tid] = thd_pdx[tid] + thd_p[tid] ;        	 	/* ptr to 1st row of d affecting s in Lx */
	      
	      for (thd_pdi2[tid] = thd_pdi1[tid] ; thd_pdi2[tid] < thd_pdend[tid] && Ls [thd_pdi2[tid]] < k2 ; (thd_pdi2[tid])++) ;
	      thd_ndrow1[tid] = thd_pdi2[tid] - thd_pdi1[tid] ;      	 	/* # rows in first part of d */
	      thd_ndrow2[tid] = thd_pdend[tid] - thd_pdi1[tid] ;     	 	/* # rows in remaining d */
	      
	      /* construct the update matrix C for this supernode d */
	      thd_ndrow3[tid] = thd_ndrow2[tid] - thd_ndrow1[tid] ;  	 	/* number of rows of C2 */
	      thd_C[tid] = ptrC;
	      ptrC = ptrC + thd_ndrow2[tid]*thd_ndrow1[tid];
	      Cused += thd_ndrow2[tid]*thd_ndrow1[tid];

	      /* increment the descendant pointer */
	      if ( tid > 0 ){
		if ( Cused < L->maxcsize) {
		  d = dsmall;
		  idescendant++;
		  dsmall = Previous[dsmall];
		  nlocalt++;
		}
		else {
		  thd_d[tid] = -1;
		}
	      }

	    }
	    else {
	      thd_d[tid] = -1;
	    }

	  }


#pragma omp parallel num_threads(THREADS_INNER) 
	  {
	    
	    /* obtain thread id */
	    int tid = 0;
	    tid = omp_get_thread_num();
	    
	    if ( thd_d[tid] > -1 ) {
	    
#ifdef REAL
	    /* dsyrk */
	    BLAS_dsyrk ("L", "N",
			thd_ndrow1[tid], thd_ndcol[tid],              		/* N, K: L1 is ndrow1-by-ndcol*/
			one,                        		                /* ALPHA:  1 */
			Lx + L_ENTRY*thd_pdx1[tid], thd_ndrow[tid],   		/* A, LDA: L1, ndrow */
			zero,                       		                /* BETA:   0 */
			thd_C[tid], thd_ndrow2[tid]) ;                		/* C, LDC: C1 */
#else
	    BLAS_zherk ("L", "N",
			thd_ndrow1[tid], thd_ndcol[tid],                        /* N, K: L1 is ndrow1-by-ndcol*/
			one,                                                    /* ALPHA:  1 */
			Lx + L_ENTRY*thd_pdx1[tid], thd_ndrow[tid],             /* A, LDA: L1, ndrow */
			zero,                                                   /* BETA:   0 */
			thd_C[tid], thd_ndrow2[tid]) ;                          /* C, LDC: C1 */
#endif
	    
	    
	    /* dgemm */
	    if (thd_ndrow3[tid] > 0)
	      {
#ifdef REAL
		BLAS_dgemm ("N", "C",
			    thd_ndrow3[tid], thd_ndrow1[tid], thd_ndcol[tid],          	/* M, N, K */
			    one,                            	                        /* ALPHA:  1 */
			    Lx + L_ENTRY*(thd_pdx1[tid] + thd_ndrow1[tid]),   	        /* A, LDA: L2 */
			    thd_ndrow[tid],                          	                /* ndrow */
			    Lx + L_ENTRY*thd_pdx1[tid],              	                /* B, LDB: L1 */
			    thd_ndrow[tid],                          	                /* ndrow */
			    zero,                           	                        /* BETA:   0 */
			    thd_C[tid] + L_ENTRY*thd_ndrow1[tid],             	        /* C, LDC: C2 */
			    thd_ndrow2[tid]) ;
#else
		BLAS_zgemm ("N", "C",
			    thd_ndrow3[tid], thd_ndrow1[tid], thd_ndcol[tid],           /* M, N, K */
			    one,                                                        /* ALPHA:  1 */
			    Lx + L_ENTRY*(thd_pdx1[tid] + thd_ndrow1[tid]),             /* A, LDA: L2 */
			    thd_ndrow[tid],                                             /* ndrow */
			    Lx + L_ENTRY*thd_pdx1[tid],                                 /* B, LDB: L1, ndrow */
			    thd_ndrow[tid],
			    zero,                                                       /* BETA:   0 */
			    thd_C[tid] + L_ENTRY*thd_ndrow1[tid],                       /* C, LDC: C2 */
			    thd_ndrow2[tid]) ;
#endif
		
	      }
	    }
	  }
	  
	  {
	  nlocalt = THREADS;

	  for ( tid=0; tid<nlocalt; tid++ ) {

	    if ( thd_d[tid] > -1 ) {

	      /* construct relative map to assemble d into s */
#pragma omp parallel for num_threads(THREADS_INNER) if ( thd_ndrow2[tid] > 64 )
	      for (i = 0 ; i < thd_ndrow2[tid] ; i++)
		RelativeMap [i] = Map [Ls [thd_pdi1[tid] + i]] ;
	      
	  
	      /* assemble C into supernode s using the relative map */
#pragma omp parallel for private ( j, i, px, q ) num_threads(8/*THREADS_INNER*/) if (thd_ndrow1[tid] > 64 )
	      for (j = 0 ; j < thd_ndrow1[tid] ; j++) {              	                /* cols k1:k2-1 */
		
		px = psx + RelativeMap [j] * nsrow ;
		for (i = j ; i < thd_ndrow2[tid] ; i++)          	                /* rows k1:n-1 */
		  {
		    q = px + RelativeMap [i] ;
		    L_ASSEMBLESUB (Lx,q, thd_C[tid], i+thd_ndrow2[tid]*j) ;		/* Lx [px + RelativeMap [i]] -= C [i + pj] ; */
		  }
	      }

	      /* prepare this supernode d for its next ancestor */
	      dnext = Next [thd_d[tid]] ;  
	      if (!repeat_supernode) {
		Lpos [thd_d[tid]] = thd_pdi2[tid] - thd_pdi[tid] ;
		if (Lpos [thd_d[tid]] < thd_ndrow[tid]) {
		  dancestor = SuperMap [Ls [thd_pdi2[tid]]] ;
		  /* place d in the link list of its next ancestor */
		  Next [thd_d[tid]] = Head [dancestor] ;
		  Head [dancestor] = thd_d[tid] ;
		}
	      }	
	    }
	  }

      }  /* end of pragma omp critical */
	  
	}
	else  {

	  supernodeUsedGPU = 1;   				/* GPU was used for this supernode*/
	  Common->ibuffer[gpuid]++;
	  Common->ibuffer[gpuid] = Common->ibuffer[gpuid]%(CHOLMOD_HOST_SUPERNODE_BUFFERS*CHOLMOD_DEVICE_STREAMS);

	  /* prepare this supernode d for its next ancestor */
	  dnext = Next [d] ;  
	  if (!repeat_supernode) {
	    Lpos [d] = pdi2 - pdi ;
	    if (Lpos [d] < ndrow) {
	      dancestor = SuperMap [Ls [pdi2]] ;
	      /* place d in the link list of its next ancestor */
	      Next [d] = Head [dancestor] ;
	      Head [dancestor] = d ;
	    }
	  }
	  
	}
	
	
      } /* end loop over descendants */
      
      

      
      /*
       *  Final Supernode Assembly
       *
       *  Sum CPU and GPU assembly's of supernode:
       *
       */
      iHostBuff = (Common->ibuffer[gpuid])%CHOLMOD_HOST_SUPERNODE_BUFFERS;
      iDevBuff = (Common->ibuffer[gpuid])%CHOLMOD_DEVICE_STREAMS;
      TEMPLATE2 ( CHOLMOD (gpu_final_assembly_gtc ))( Common, gpu_p, Lx, &iHostBuff, &iDevBuff, psx, nscol, nsrow, supernodeUsedGPU, gpuid );
     

  
      /*
       *  Cholesky Factorization
       *
       *  Factorize diagonal block of spuernode s in LL' in the following steps:
       *  1. perform dpotrf
       *
       */
      nscol2 = (repeat_supernode) ? (nscol3) : (nscol) ;
      if ( !(supernodeUsedGPU) || !TEMPLATE2 (CHOLMOD (gpu_lower_potrf_gtc))(Common, gpu_p, Lx, &info, nscol2, nsrow, psx, gpuid))
      {

        supernodeUsedGPU = 0;

#ifdef REAL
        LAPACK_dpotrf ("L",
                       nscol2,                    	/* N: nscol2 */
                       Lx + L_ENTRY*psx, nsrow,    	/* A, LDA: S1, nsrow */
                       info) ;                    	/* INFO */    
#else
        LAPACK_zpotrf ("L",
        	       nscol2,                     /* N: nscol2 */
                       Lx + L_ENTRY*psx, nsrow,    /* A, LDA: S1, nsrow */
                       info) ;                     /* INFO */
#endif


      }




      /* check if the matrix is not positive definite */
      if (repeat_supernode)
      {
        /* the leading part has been refactorized; it must have succeeded */
        info = 0 ;
  
        /* zero out the rest of this supernode */
        p = psx + nsrow * nscol3 ;
        pend = psx + nsrow * nscol ;            
        for ( ; p < pend ; p++)
        {
          L_CLEAR (Lx,p) ;				/* Lx [p] = 0 ; */
        }
      }




      /* info is set to one in LAPACK_*potrf if blas_ok is FALSE.  It is
       * set to zero in dpotrf/zpotrf if the factorization was successful. */
      if (CHECK_BLAS_INT && !Common->blas_ok)
      {
        ERROR (CHOLMOD_TOO_LARGE, "problem too large for the BLAS") ;
      }




      /* check if the matrix is not positive definite */
      if (info != 0)
      {

        /* Matrix is not positive definite.  dpotrf/zpotrf do NOT report an
         * error if the diagonal of L has NaN's, only if it has a zero. 
         */
        if (Common->status == CHOLMOD_OK)
        {
          ERROR (CHOLMOD_NOT_POSDEF, "matrix not positive definite") ;
        }


        /* L->minor is the column of L that contains a zero or negative
         * diagonal term. 
         */
        L->minor = k1 + info - 1 ;


        /* clear the link lists of all subsequent supernodes */
        for (ss = s+1 ; ss < nsuper ; ss++)
        {
          Head [ss] = EMPTY ;
        }


        /* zero this supernode, and all remaining supernodes */
        pend = L->xsize ;
        for (p = psx ; p < pend ; p++)
        {
          /* Lx [p] = 0. ; */
          L_CLEAR (Lx,p) ;
        }


        /* If L is indefinite, it still contains useful information.
         * Supernodes 0 to s-1 are valid, similar to MATLAB [R,p]=chol(A),
         * where the 1-based p is identical to the 0-based L->minor.  Since
         * L->minor is in the current supernode s, it and any columns to the
         * left of it in supernode s are also all zero.  This differs from
         * [R,p]=chol(A), which contains nonzero rows 1 to p-1.  Fix this
         * by setting repeat_supernode to TRUE, and repeating supernode s.
         *
         * If Common->quick_return_if_not_posdef is true, then the entire
         * supernode s is not factorized; it is left as all zero.
         */
         if (info == 1 || Common->quick_return_if_not_posdef)
         {
           /* If the first column of supernode s contains a zero or
            * negative diagonal entry, then it is already properly set to
            * zero.  Also, info will be 1 if integer overflow occured in
            * the BLAS. */
            Head [s] = EMPTY ;

            /* finalize GPU */
            CHOLMOD (gpu_end) (Common);
            /*return (Common->status >= CHOLMOD_OK) ;*/

          }
          else
          {
            /* Repeat supernode s, but only factorize it up to but not
             * including the column containing the problematic diagonal
             * entry. */
             repeat_supernode = TRUE ;
             s-- ;
             nscol3 = info - 1 ;
          }
    
      } /* end if info */



      /*
       *  Triangular Solve
       *
       *  Compute the subdiagonal block in the following steps:
       *  1. perform dtrsm
       *  2. copy result back into factor Lx 
       *  3. prepare next supernode
       *
       */
      nsrow2 = nsrow - nscol2 ;
      if (nsrow2 > 0)
      { 
        /* The current supernode is columns k1 to k2-1 of L.  Let L1 be the
         * diagonal block (factorized by dpotrf/zpotrf above; rows/cols
         * k1:k2-1), and L2 be rows k2:n-1 and columns k1:k2-1 of L.  The
         * triangular system to solve is L2*L1' = S2, where S2 is
         * overwritten with L2.  More precisely, L2 = S2 / L1' in MATLAB
         * notation.
         */
         if ( !(supernodeUsedGPU) || !TEMPLATE2 (CHOLMOD(gpu_triangular_solve_gtc)) (Common, gpu_p, Lx, nsrow2, nscol2, nsrow, psx ,gpuid) )
         {
#ifdef REAL
           BLAS_dtrsm ("R", "L", "C", "N",
                       nsrow2, nscol2,                 	/* M, N */
                       one,                            	/* ALPHA: 1 */
                       Lx + L_ENTRY*psx, nsrow,        	/* A, LDA: L1, nsrow */
                       Lx + L_ENTRY*(psx + nscol2),    	/* B, LDB, L2, nsrow */
                       nsrow) ;
#else
           BLAS_ztrsm ("R", "L", "C", "N",
                       nsrow2, nscol2,                 /* M, N */
                       one,                            /* ALPHA: 1 */
                       Lx + L_ENTRY*psx, nsrow,        /* A, LDA: L1, nsrow */
                       Lx + L_ENTRY*(psx + nscol2),    /* B, LDB, L2, nsrow */
                       nsrow) ;
#endif

    
         }


         if (CHECK_BLAS_INT && !Common->blas_ok)
         {
           ERROR (CHOLMOD_TOO_LARGE, "problem too large for the BLAS") ;
         }


         /* prepare next supernode */  
         if (!repeat_supernode)
         {
           /* Lpos [s] is offset of first row of s affecting its parent */
           Lpos [s] = nscol ;
           sparent = SuperMap [Ls [psi + nscol]] ;
           /* place s in link list of its parent */
           Next [s] = Head [sparent] ;
           Head [sparent] = s ;
         }
      }
      /* copy supernode back to factor Lx anyways */
      else
      {
        TEMPLATE2 ( CHOLMOD (gpu_copy_supernode_gtc) )( Common, gpu_p, Lx, psx, nscol, nscol2, nsrow, supernodeUsedGPU, iHostBuff, gpuid);
      }

      Head [s] = EMPTY ;  /* link list for supernode s no longer needed */


    } /* end loop over supenodes */


  } /* end loop over levels */



}


/*
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
*/











