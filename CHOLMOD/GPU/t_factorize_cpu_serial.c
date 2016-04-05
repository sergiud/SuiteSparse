/* ========================================================================== */
/* === GPU/t_factorize_cpu_serial.c ========================================= */
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
 *   t_factorize_cpu_serial
 *
 * Description:
 *   Contains functions for factorization
 *   of the CPU algorithm. (serial version)
 *
 */


/* includes */
#include <string.h>
#include <time.h>
#include "cholmod_template.h"

#ifdef MKLROOT
#include "mkl.h"
#endif








/*
 * Function:
 *   gpu_factorize_cpu_serial
 *
 * Description:
 *   Factorizes entire elimination tree on the CPU (serial version)
 *   Returns 0 if matrix not positive-definite, 1 otherwise.
 *
 */
int TEMPLATE2 (CHOLMOD (gpu_factorize_cpu_serial))
  (
    cholmod_common *Common,
    cholmod_factor *L,
    cholmod_global_pointers *gb_p,
    cholmod_cpu_pointers *cpu_p,
    cholmod_tree_pointers *tree_p,
    cholmod_profile_pointers *prof_p,
    int deviceid
   )
{
  /* local variables */
  int i, j, k, numThreads;
  Int px, pk, pf, p, q, d, s, n, ss, nsuper, ndrow, ndrow1, ndrow2, ndrow3, ndcol, nsrow, nsrow2, nscol, nscol2, nscol3, 
      kd1, kd2, k1, k2, psx, psi, pdx, pdx1, pdi, pdi1, pdi2, pdend, psend, pfend, pend, dancestor, sparent, imap, 
      start, end, super, dnext, info, repeat_supernode, Apacked, Fpacked, stype;
  Int *Ls, *Lpi, *Lpx, *Lpos, *Fp, *Fi, *Fnz, *Ap, *Ai, *Anz, *Super, *Map, *RelativeMap, *SuperMap, *Head, *Next, *Next_save, *Lpos_save;
  double *Lx, *Ax, *Az, *Fx, *Fz, *C, *beta, *tstart, *tend, *syrk_time, *gemm_time, *potrf_time, *trsm_time, 
	  *syrk_flops, *gemm_flops, *potrf_flops, *trsm_flops;
  double one[2] = {1.0, 0.0}, zero[2] = {0.0, 0.0};
  double tstart1;





  /*
   * Set variables & pointers
   */
  /* set host variables */
  n 			= L->n;
  numThreads		= Common->ompNumThreads;
  nsuper		= L->nsuper;
  repeat_supernode 	= FALSE ;
  Apacked       	= cpu_p->Apacked;
  Fpacked       	= cpu_p->Fpacked;
  stype         	= cpu_p->stype;
  beta          	= cpu_p->beta;

  /* set host pointers */
  Ls            = cpu_p->Ls;
  Lpi           = cpu_p->Lpi;
  Lpx 		= L->px;
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
  Lx            = cpu_p->Lx;
  Ax            = cpu_p->Ax;
  Az            = cpu_p->Az;
  Fx            = cpu_p->Fx;
  Fz            = cpu_p->Fz;
  C             = cpu_p->C;

  /* set timer pointers */
  tstart        = prof_p->f_start[deviceid];
  tend          = prof_p->f_end[deviceid];
  syrk_time     = prof_p->syrk_time[deviceid];
  gemm_time     = prof_p->gemm_time[deviceid];
  potrf_time    = prof_p->potrf_time[deviceid];
  trsm_time     = prof_p->trsm_time[deviceid];
  syrk_flops    = prof_p->syrk_flop[deviceid];
  gemm_flops    = prof_p->gemm_flop[deviceid];
  potrf_flops   = prof_p->potrf_flop[deviceid];
  trsm_flops    = prof_p->trsm_flop[deviceid];




#ifdef MKLROOT
  /* set mkl threads */  
  mkl_set_num_threads(numThreads);
#else
  openblas_set_num_threads(numThreads);
#endif




  /* clear global flops */
  CLEAR1(syrk_flops,0);
  CLEAR1(gemm_flops,0);
  CLEAR1(potrf_flops,0);
  CLEAR1(trsm_flops,0);

 
    TIMER_START(tstart,0);     


    /* loop over supernodes */
    for(super = 0; super < nsuper; super++)
    {


      /* get supernode dimensiosn */
      s = super;
      k1 = Super [s] ;            		/* s contains columns k1 to k2-1 of L */
      k2 = Super [s+1] ;
      nscol = k2 - k1 ;           		/* # of columns in all of s */
      psi = Lpi [s] ;             		/* pointer to first row of s in Ls */
      psx = Lpx [s] ;             		/* pointer to first row of s in Lx */
      psend = Lpi [s+1] ;         		/* pointer just past last row of s in Ls */
      nsrow = psend - psi ;       		/* # of rows in all of s */
      pend = psx + nsrow * nscol ;      	/* s is nsrow-by-nscol */
      pk = psx ;
      SUM(potrf_flops,0,(double)(nscol*nscol*nscol/3.0));
      SUM(trsm_flops,0,(double)((nsrow-nscol)*nscol*nscol));

  
  


      TIMER_START(tstart,3);
      /* construct the scattered Map for supernode s */
      #pragma omp parallel for num_threads(numThreads) if ( nsrow > 128 )
      for (k = 0 ; k < nsrow ; k++)
        Map [Ls [psi + k]] = k ;
  



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
                L_ASSIGN (Lx,(imap+(psx+(k-k1)*nsrow)), Ax,Az,p) 		/* Lx [Map [i] + pk] = Ax [p] ; */;
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
          L_ASSEMBLE (Lx,pk, beta) ;		/* Lx [pk] += beta [0] ; */
          pk += nsrow + 1 ;       		/* advance to the next diagonal entry */
        }
      }

      TIMER_END(tstart,tend,3);      



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




      dnext = Head[s];
      /* loop over descendant d of supernode s  */
      while( (dnext != EMPTY)  || (0) ) 
      {

        d = dnext;

  
        /* get the size of supernode d */
        kd1 = Super [d] ;       		/* d contains cols kd1 to kd2-1 of L */
        kd2 = Super [d+1] ;
        ndcol = kd2 - kd1 ;     		/* # of columns in all of d */
        pdi = Lpi [d] ;         		/* pointer to first row of d in Ls */
        pdx = Lpx [d] ;         		/* pointer to first row of d in Lx */
        pdend = Lpi [d+1] ;     		/* pointer just past last row of d in Ls */
        ndrow = pdend - pdi ;   		/* # rows in all of d */
  
        /* find the range of rows of d that affect rows k1 to k2-1 of s */
        p = Lpos [d] ;          		/* offset of 1st row of d affecting s */
        pdi1 = pdi + p ;        		/* ptr to 1st row of d affecting s in Ls */
        pdx1 = pdx + p ;   		        /* ptr to 1st row of d affecting s in Lx */
  
        for (pdi2 = pdi1 ; pdi2 < pdend && Ls [pdi2] < k2 ; (pdi2)++) ;
        ndrow1 = pdi2 - pdi1 ;      		/* # rows in first part of d */
        ndrow2 = pdend - pdi1 ;     		/* # rows in remaining d */
  
        /* construct the update matrix C for this supernode d */
        ndrow3 = ndrow2 - ndrow1 ;  /* number of rows of C2 */
        SUM(syrk_flops,0,(double)(ndrow1*ndrow1*ndcol));
        SUM(gemm_flops,0,(double)(2*(ndrow2-ndrow1)*ndrow1*ndcol));



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

        TIMER_START1(tstart1);          
#ifdef REAL
        /* dsyrk */
        BLAS_dsyrk ("L", "N",
        	    ndrow1, ndcol,              /* N, K: L1 is ndrow1-by-ndcol*/
        	    one,                        /* ALPHA:  1 */
        	    Lx + L_ENTRY*pdx1, ndrow,   /* A, LDA: L1, ndrow */
        	    zero,                       /* BETA:   0 */
        	    C, ndrow2) ;                /* C, LDC: C1 */
#else
        BLAS_zherk ("L", "N",
        	    ndrow1, ndcol,              /* N, K: L1 is ndrow1-by-ndcol*/
        	    one,                        /* ALPHA:  1 */
        	    Lx + L_ENTRY*pdx1, ndrow,   /* A, LDA: L1, ndrow */
        	    zero,                       /* BETA:   0 */
        	    C, ndrow2) ;                /* C, LDC: C1 */
#endif
        TIMER_END1(tstart1,syrk_time,0);





        /* dgemm */
        TIMER_START1(tstart1);
        if (ndrow3 > 0)
        {
#ifdef REAL
          BLAS_dgemm ("N", "C",
          	      ndrow3, ndrow1, ndcol,          /* M, N, K */
          	      one,                            /* ALPHA:  1 */
          	      Lx + L_ENTRY*(pdx1 + ndrow1),   /* A, LDA: L2 */
          	      ndrow,                          /* ndrow */
          	      Lx + L_ENTRY*pdx1,              /* B, LDB: L1 */
          	      ndrow,                          /* ndrow */
          	      zero,                           /* BETA:   0 */
          	      C + L_ENTRY*ndrow1,             /* C, LDC: C2 */
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
        TIMER_END1(tstart1,gemm_time,0);





        TIMER_START(tstart,4);
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
            L_ASSEMBLESUB (Lx,q, C, i+ndrow2*j) ;	/* Lx [px + RelativeMap [i]] -= C [i + pj] ; */
          }
        }
        TIMER_END(tstart,tend,4);




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


      } /* end of descendant supernode loop */




      /*
       *  Cholesky Factorization
       *
       *  Factorize diagonal block of spuernode s in LL' in the following steps:
       *  1. perform dpotrf
       *
       */
      TIMER_START1(tstart1);
      nscol2 = (repeat_supernode) ? (nscol3) : (nscol) ;
#ifdef REAL
      LAPACK_dpotrf ("L",
                     nscol2,                    	/* N: nscol2 */
                     Lx + L_ENTRY*psx, nsrow,    	/* A, LDA: S1, nsrow */
                     info) ;                    	/* INFO */
#else
      LAPACK_zpotrf ("L",
                     nscol2,                     	/* N: nscol2 */
                     Lx + L_ENTRY*psx, nsrow,    	/* A, LDA: S1, nsrow */
                     info) ;                     	/* INFO */
#endif
        TIMER_END1(tstart1,potrf_time,0);

    


      /* check if the matrix is not positive definite */
      if (repeat_supernode)
      {
        /* the leading part has been refactorized; it must have succeeded */
        info = 0 ;

        /* zero out the rest of this supernode */
        p = psx + nsrow * nscol3 ;
        pend = psx + nsrow * nscol ;            /* s is nsrow-by-nscol */
        for ( ; p < pend ; p++)
        {
          /* Lx [p] = 0 ; */
          L_CLEAR (Lx,p) ;
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
         * error if the diagonal of L has NaN's, only if it has a zero. */
        if (Common->status == CHOLMOD_OK)
        {
          ERROR (CHOLMOD_NOT_POSDEF, "matrix not positive definite") ;
        }

        /* L->minor is the column of L that contains a zero or negative
         * diagonal term. */
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
            return Common->status;/*(Common->status >= CHOLMOD_OK)*/
  
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
    
      }




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

         /* dtrsm */
	 TIMER_START1(tstart1);
#ifdef REAL
         BLAS_dtrsm ("R", "L", "C", "N",
                     nsrow2, nscol2,                 	/* M, N */
                     one,                            	/* ALPHA: 1 */
                     Lx + L_ENTRY*psx, nsrow,        	/* A, LDA: L1, nsrow */
                     Lx + L_ENTRY*(psx + nscol2),    	/* B, LDB, L2, nsrow */
                     nsrow) ;
#else
         BLAS_ztrsm ("R", "L", "C", "N",
                     nsrow2, nscol2,                 	/* M, N */
                     one,                            	/* ALPHA: 1 */
                     Lx + L_ENTRY*psx, nsrow,        	/* A, LDA: L1, nsrow */
                     Lx + L_ENTRY*(psx + nscol2),    	/* B, LDB, L2, nsrow */
                     nsrow) ;
#endif
        TIMER_END1(tstart1,trsm_time,0);

  


         if (CHECK_BLAS_INT && !Common->blas_ok)
         {
           ERROR (CHOLMOD_TOO_LARGE, "problem too large for the BLAS") ;
         }
  

         /* prepare for next supernode */
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


      Head [s] = EMPTY ;  /* link list for supernode s no longer needed */



        if (repeat_supernode)
        {
          /* matrix is not positive definite; finished clean-up for supernode
           * containing negative diagonal */
          return Common->status;/*(Common->status >= CHOLMOD_OK)*/
        }


    } /* end loop over supenodes */
    TIMER_END(tstart,tend,0);    





  /* print overall benchmarks */
  PRINTF("\n\nElimination tree benchmarks:\n");
  PRINTF("\n- time -\n");
  PRINTFV("total:        %f\n",tend[0]);
  PRINTFV("initLx:       %f\n",tend[3]);
  PRINTFV("dsyrk:        %f\n",syrk_time[0]);
  PRINTFV("dgemm:        %f\n",gemm_time[0]);
  PRINTFV("assembly:     %f\n",tend[4]);
  PRINTFV("dpotrf:       %f\n",potrf_time[0]);
  PRINTFV("dtrsm:        %f\n",trsm_time[0]);
  PRINTF("\n- flops -\n");
  PRINTFV("dsyrk:        %f\n",1.0e-9*syrk_flops[0]/syrk_time[0]);
  PRINTFV("dgemm:        %f\n",1.0e-9*gemm_flops[0]/gemm_time[0]);
  PRINTFV("dpotrf:       %f\n",1.0e-9*potrf_flops[0]/potrf_time[0]);
  PRINTFV("dtrsm:        %f\n",1.0e-9*trsm_flops[0]/trsm_time[0]);
  PRINTF("\n");

  

  /* return ok */
  return Common->status;/*(Common->status >= CHOLMOD_OK)*/
}


/*
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
*/
/*
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
*/











