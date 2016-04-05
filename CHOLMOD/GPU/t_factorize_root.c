/* ========================================================================== */
/* === GPU/t_factorize_root.c =============================================== */
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
 *   t_factorize_root
 *
 * Description:
 *   Contains functions for factorization of the root algorithm.
 *   Returns 1 if matrix not positive-definite, 0 otherwise.
 *
 */


/* includes */
#include <string.h>
#include <time.h>


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
#define L_ASSEMBLESUB(Lx,q,C,p)         Lx [2*(q)  ] -= C [2*(p)  ] ; \
                                        Lx [2*(q)+1] -= C [2*(p)+1] ;

#ifdef COMPLEX
/* A, F, L, and C are all complex */
#define L_ASSIGN(Lx,q, Ax,Az,p)         Lx [2*(q)  ] = Ax [2*(p)  ] ; \
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
 *   gpu_factorize_root
 *
 * Description:
 *   Factorizes top-of-tree of elimination tree, where
 *   the subtree does not fit the GPU. Utilizes a hybrid algorithm
 *   presented at GTC14.
 *   Returns 0 if matrix not positive definite, 1 otherwise.
 *
 */
int TEMPLATE2 (CHOLMOD (gpu_factorize_root))
  (
    cholmod_common *Common,
    cholmod_factor *L,
    cholmod_gpu_pointers *gpu_p,
    cholmod_cpu_pointers *cpu_p,
    cholmod_tree_pointers *tree_p,
    Int subtree
   )
{
#ifdef SUITESPARSE_CUDA
  /* local variables */
  int i, j, k, gpuid, numThreads;
  Int px, pk, pf, p, q, d, s, n, ss, ndrow, ndrow1, ndrow2, ndrow3, ndcol, nsrow, nsrow2, nscol, nscol2, nscol3, 
      kd1, kd2, k1, k2, psx, psi, pdx, pdx1, pdi, pdi1, pdi2, pdend, psend, pfend, pend, dancestor, sparent, imap, 
      start, end, node, level, idescendant, dnext, dlarge, iHostBuff, iDevBuff, skips, dsmall, tail, info,
      GPUavailable, mapCreatedOnGpu, supernodeUsedGPU, repeat_supernode, Apacked, Fpacked, stype;
  Int *Ls, *Lpi, *Lpx, *Lpos, *Fp, *Fi, *Fnz, *Ap, *Ai, *Anz, *Super, *Map, *RelativeMap, *SuperMap, *Head, *Next, *Next_save, *Previous, *Lpos_save,
      *ndescendants, *supernode_levels, *supernode_levels_ptrs, *supernode_levels_subtree_ptrs, *supernode_num_levels;
  double *Lx, *Ax, *Az, *Fx, *Fz, *C, *beta;
  double one[2] = {1.0, 0.0}, zero[2] = {0.0, 0.0};
  cudaError_t cuErr;





  /*
   * Set variables & pointers
   */
  /* set host variables */
  n	 		= L->n;
  nscol3 		= 0;
  numThreads		= Common->ompNumThreads;
  gpu_p->gpuid 		= 0;
  gpuid 		= gpu_p->gpuid;
  repeat_supernode 	= FALSE;
  Apacked               = cpu_p->Apacked;
  Fpacked               = cpu_p->Fpacked;
  stype                 = cpu_p->stype;
  beta                  = cpu_p->beta;

  /* set host pointers */
  Ls            = cpu_p->Ls;
  Lpi           = cpu_p->Lpi;
  Lpx           = L->px;
  Lpos          = cpu_p->Lpos;
  Fp            = cpu_p->Fp;
  Fi            = cpu_p->Fi;
  Fnz           = cpu_p->Fnz;
  Ap            = cpu_p->Ap;
  Ai            = cpu_p->Ai;
  Anz           = cpu_p->Anz;
  Super         = cpu_p->Super;
  Map           = cpu_p->Map;
  RelativeMap   = cpu_p->RelativeMap;
  SuperMap      = cpu_p->SuperMap;
  Head          = cpu_p->Head;
  Next          = cpu_p->Next;
  Next_save     = cpu_p->Next_save;
  Lpos_save     = cpu_p->Lpos_save;
  Previous	= cpu_p->Previous;
  Lx            = cpu_p->Lx;
  Ax            = cpu_p->Ax;
  Az            = cpu_p->Az;
  Fx            = cpu_p->Fx;
  Fz            = cpu_p->Fz;
  C             = cpu_p->C;

  /* set tree pointers */
  ndescendants              	= tree_p->ndescendants;
  supernode_levels              = tree_p->supernode_levels;
  supernode_levels_ptrs         = tree_p->supernode_levels_ptrs;
  supernode_levels_subtree_ptrs  = tree_p->supernode_levels_subtree_ptrs;
  supernode_num_levels          = tree_p->supernode_num_levels;





  /* initialize GPU */
  TEMPLATE2 (CHOLMOD (gpu_init_root))(Common, gpu_p, L, Lpi, L->nsuper, n,gpuid);




 
  /* loop over levels in subtree */
  for(level = 0; level < supernode_num_levels[subtree]; level++)
  {  

    start = supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+level];
    end = supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+level+1];

    /* loop over supernodes */
    for(node = start; node < end; node++)
    {


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



  
      /* construct the scattered Map for supernode s */
      #pragma omp parallel for num_threads(numThreads) if ( nsrow > 128 )
      for (k = 0 ; k < nsrow ; k++)
        Map [Ls [psi + k]] = k ;
  



      /* reorder descendants in supernode by descreasing size */
      TEMPLATE2 (CHOLMOD (gpu_reorder_descendants_root))(Common, gpu_p, Lpi, Lpos, Super, Head, &tail, Next, Previous,
                                                        &ndescendants[0], &mapCreatedOnGpu, s, gpuid );




      /* copy matrix into supernode s (lower triangular part only) */
      #pragma omp parallel for private ( p, pend, pfend, pf, i, j, imap, q ) num_threads(numThreads) if ( k2-k1 > 64 )
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
      skips = 0;




      /* loop over descendants d of supernode s */
      while( (idescendant < ndescendants[0]) )
      {
  
        iHostBuff = (Common->ibuffer[gpuid]) % CHOLMOD_HOST_SUPERNODE_BUFFERS;


        /* get next descendant */  
        if ( idescendant > 0 )  {
          if ( GPUavailable == -1 || skips > 0) {
            d = dsmall;
            dsmall = Previous[dsmall];
            (skips)--;
          }
          else {
            cuErr = cudaEventQuery( Common->updateCBuffersFree[gpuid][iHostBuff] );
            if ( cuErr == cudaSuccess ) {
              d = dlarge;
              dlarge = Next[dlarge];
              GPUavailable = 1;
              skips = 0;
            }
            else {
              d = dsmall;
              dsmall = Previous[dsmall];
              GPUavailable = 0;
              skips = CHOLMOD_GPU_SKIP;
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
                TEMPLATE2 ( CHOLMOD (gpu_initialize_supernode_root))( Common, gpu_p, nscol, nsrow, psi, gpuid );
                mapCreatedOnGpu = 1;
              }
            }
          else {
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
        if ( GPUavailable!=1 || !TEMPLATE2 (CHOLMOD (gpu_updateC_root)) (Common, gpu_p, Lx, C, ndrow1, ndrow2, ndrow, ndcol, nsrow, pdx1, pdi1, gpuid) )
        {

#ifdef REAL
          /* dsyrk */
          BLAS_dsyrk ("L", "N",
          	      ndrow1, ndcol,              		/* N, K: L1 is ndrow1-by-ndcol*/
          	      one,                        		/* ALPHA:  1 */
          	      Lx + L_ENTRY*pdx1, ndrow,   		/* A, LDA: L1, ndrow */
          	      zero,                       		/* BETA:   0 */
          	      C, ndrow2) ;                		/* C, LDC: C1 */
#else
          BLAS_zherk ("L", "N",
                      ndrow1, ndcol,              /* N, K: L1 is ndrow1-by-ndcol*/
                      one,                        /* ALPHA:  1 */
                      Lx + L_ENTRY*pdx1, ndrow,   /* A, LDA: L1, ndrow */
                      zero,                       /* BETA:   0 */
                      C, ndrow2) ;                /* C, LDC: C1 */
#endif


          /* dgemm */
          if (ndrow3 > 0)
          {
#ifdef REAL
            BLAS_dgemm ("N","C",
            	  	ndrow3, ndrow1, ndcol,          	/* M, N, K */
            		one,                            	/* ALPHA:  1 */
            		Lx + L_ENTRY*(pdx1 + ndrow1),   	/* A, LDA: L2 */
            		ndrow,                          	/* ndrow */
            		Lx + L_ENTRY*pdx1,              	/* B, LDB: L1 */
            		ndrow,                          	/* ndrow */
            		zero,                           	/* BETA:   0 */
            		C + L_ENTRY*ndrow1,             	/* C, LDC: C2 */
            		ndrow2) ;
#else
            BLAS_zgemm ("N", "C",
                        ndrow3, ndrow1, ndcol,          /* M, N, K */
                        one,                            /* ALPHA:  1 */
                        Lx + L_ENTRY*(pdx1 + ndrow1),   /* A, LDA: L2 */
                        ndrow,                          /* ndrow */
                        Lx + L_ENTRY*pdx1,              /* B, LDB: L1, ndrow */
                        ndrow,
                        zero,                           /* BETA:   0 */
                        C + L_ENTRY*ndrow1,             /* C, LDC: C2 */
                        ndrow2) ;
#endif

          }


          /* construct relative map to assemble d into s */
          #pragma omp parallel for num_threads(numThreads) if ( ndrow2 > 64 )
            for (i = 0 ; i < ndrow2 ; i++)
              RelativeMap [i] = Map [Ls [pdi1 + i]] ;


          /* assemble C into supernode s using the relative map */
          #pragma omp parallel for private ( j, i, px, q ) num_threads(numThreads) if (ndrow1 > 64 )
              for (j = 0 ; j < ndrow1 ; j++)              	/* cols k1:k2-1 */
              {
                px = psx + RelativeMap [j] * nsrow ;
                for (i = j ; i < ndrow2 ; i++)          	/* rows k1:n-1 */
                {
                  q = px + RelativeMap [i] ;
                  L_ASSEMBLESUB (Lx,q, C, i+ndrow2*j) ;		/* Lx [px + RelativeMap [i]] -= C [i + pj] ; */
                }
              }

          }
          else
          {
            supernodeUsedGPU = 1;   				/* GPU was used for this supernode*/
            Common->ibuffer[gpuid]++;
            Common->ibuffer[gpuid] = Common->ibuffer[gpuid]%(CHOLMOD_HOST_SUPERNODE_BUFFERS*CHOLMOD_DEVICE_STREAMS);
          }




        /* prepare this supernode d for its next ancestor */
        dnext = Next [d] ;  
        if (!repeat_supernode)
        {
          Lpos [d] = pdi2 - pdi ;
          if (Lpos [d] < ndrow)
          {
            dancestor = SuperMap [Ls [pdi2]] ;
            /* place d in the link list of its next ancestor */
            Next [d] = Head [dancestor] ;
            Head [dancestor] = d ;
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
      TEMPLATE2 ( CHOLMOD (gpu_final_assembly_root ))( Common, gpu_p, Lx, &iHostBuff, &iDevBuff, psx, nscol, nsrow, supernodeUsedGPU, gpuid );
     


  
      /*
       *  Cholesky Factorization
       *
       *  Factorize diagonal block of spuernode s in LL' in the following steps:
       *  1. perform dpotrf
       *
       */
      nscol2 = (repeat_supernode) ? (nscol3) : (nscol) ;
      if ( !(supernodeUsedGPU) || !TEMPLATE2 (CHOLMOD (gpu_lower_potrf_root))(Common, gpu_p, Lx, &info, nscol2, nsrow, psx, gpuid))
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
        for (ss = s+1 ; ss < L->nsuper ; ss++)
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
            return Common->status;/*(Common->status >= CHOLMOD_OK) ;*/

          }
          else
          {
            /* Repeat supernode s, but only factorize it up to but not
             * including the column containing the problematic diagonal
             * entry. */
             repeat_supernode = TRUE ;
             s-- ;
             nscol3 = info - 1 ;
             continue ;
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
         if ( !(supernodeUsedGPU) || !TEMPLATE2 (CHOLMOD(gpu_triangular_solve_root)) (Common, gpu_p, Lx, nsrow2, nscol2, nsrow, psx ,gpuid) )
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
        TEMPLATE2 ( CHOLMOD (gpu_copy_supernode_root) )( Common, gpu_p, Lx, psx, nscol, nscol2, nsrow, supernodeUsedGPU, iHostBuff, gpuid);
      }

      Head [s] = EMPTY ;  /* link list for supernode s no longer needed */


      if (repeat_supernode)
      {
        /* matrix is not positive definite; finished clean-up for supernode
         * containing negative diagonal */
        return Common->status;/*(Common->status >= CHOLMOD_OK)*/
      }


    } /* end loop over supenodes */


  } /* end loop over levels */

#endif

  /* return ok */
  return Common->status;/*(Common->status >= CHOLMOD_OK) ;*/

}


/*
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
*/











