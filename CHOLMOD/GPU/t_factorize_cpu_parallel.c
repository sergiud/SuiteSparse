/* ========================================================================== */
/* === GPU/t_factorize_cpu_parallel ========================================= */
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
 *   t_factorize_cpu_parallel
 *
 * Description:
 *   Contains functions for factorization
 *   of the CPU algorithm. (parallel version)
 *
 */


/* includes */
#include <string.h>
#include <time.h>

#ifdef MKLROOT
#include "mkl.h"
#endif










/*
 * Function:
 *   gpu_factorize_cpu_parallel
 *
 * Description:
 *   Factorizes subtree of elimination tree on the CPU (parallel version)
 *   Returns 0 if matrix not positive-definite, 1 otherwise.
 *
 */
int TEMPLATE2 (CHOLMOD (gpu_factorize_cpu_parallel))
  (
    cholmod_common *Common,
    cholmod_factor *L,
    cholmod_global_pointers *gb_p,
    cholmod_cpu_pointers *cpu_p,
    cholmod_tree_pointers *tree_p,
    cholmod_profile_pointers *prof_p,
    int deviceid,
    Int subtree
   )
{
  /* local variables */
  int i, j, super_count, desc_count, syrk_count, gemm_count, potrf_count, trsm_count, nbatch, numThreads;
  Int n, level, counter, Apacked, Fpacked, stype, iinfo[gb_p->maxbatch];;
  Int *Ls, *Lpi, *Lpx, *Lpos, *Fp, *Fi, *Fnz, *Ap, *Ai, *Anz, *Super, *SuperMap, *Head, *Next,
      *h_Map, *h_C, *ndescendants, *supernode_batch, *supernode_levels, *supernode_levels_ptrs, 
      *supernode_levels_subtree_ptrs, *supernode_num_levels;
  double *Lx, *Ax, *Az, *Fx, *Fz, *beta, *tstart, *tend, *syrk_time, *gemm_time, *potrf_time, 
	 *trsm_time, *syrk_flops, *gemm_flops, *potrf_flops, *trsm_flops;
  struct cholmod_super_t super[gb_p->maxbatch];
  struct cholmod_desc_t desc[gb_p->maxndesc];
  struct cholmod_syrk_t syrk[gb_p->maxndesc];
  struct cholmod_gemm_t gemm[gb_p->maxndesc];
  struct cholmod_potrf_t potrf[gb_p->maxbatch];
  struct cholmod_trsm_t trsm[gb_p->maxbatch];
  double tstart1;






  /*
   * Set variables & pointers
   */
  /* set host variables */
  n 			= L->n;
  numThreads 		= Common->ompNumThreads;
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
  SuperMap      = cpu_p->SuperMap;
  Head          = cpu_p->Head;
  Next          = cpu_p->Next;
  Lx            = cpu_p->Lx;
  Ax            = cpu_p->Ax;
  Az            = cpu_p->Az;
  Fx            = cpu_p->Fx;
  Fz            = cpu_p->Fz;
  h_Map		= cpu_p->Map;
  h_C		= cpu_p->C;

  /* set tree pointers */
  ndescendants              	= tree_p->ndescendants;
  supernode_batch              	= tree_p->supernode_batch;
  supernode_levels              = tree_p->supernode_levels;
  supernode_levels_ptrs         = tree_p->supernode_levels_ptrs;
  supernode_levels_subtree_ptrs  = tree_p->supernode_levels_subtree_ptrs;
  supernode_num_levels          = tree_p->supernode_num_levels;

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

  /* clear global flops */
  CLEAR1(syrk_flops,0);
  CLEAR1(gemm_flops,0);
  CLEAR1(potrf_flops,0);
  CLEAR1(trsm_flops,0);
  PRINTF("\n\nlevel:\tbatch:\tnsuper:\tsyrkTime:\tgemmTime:\tpotrfTime:\ttrsmTime:\tsyrkFlops:\tgemmFlops:\tpotrfFlops:\ttrsmFlops:\n");





  /* loop over levels in subtree */
  for(level = 0; level < supernode_num_levels[subtree]; level++)
  {


    /* set batching parameters */
    Int node = 0;
    Int start = supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+level];    	/* starting supernode of level */
    Int end   = supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+level+1];  	/* ending supernode of level */
    nbatch = supernode_batch[supernode_levels_subtree_ptrs[subtree] + level];          	/* size of batch for current level */


    /* clear local flops */
    CLEAR1(syrk_flops,1);
    CLEAR1(gemm_flops,1);
    CLEAR1(potrf_flops,1);
    CLEAR1(trsm_flops,1);


    /* loop over groups of batches */
    while(node < (end - start)) {


      /* check if batch overexceeds */
      if( (node + nbatch) > (end-start) )
        nbatch = (end-start) - node;


      TIMER_START(tstart,0);      

      /* clear counters */
      super_count = 0;
      desc_count  = 0;
      syrk_count  = 0;
      gemm_count  = 0;
      potrf_count = 0;
      trsm_count  = 0;
      counter     = 0;





      /*
       *  Store supernode dimensions
       *
       *  Stores dimensions of supernodes in the following order:
       *    1. potrf's to be sent to cuBlas
       *    2. trsm's to be sent to cuBlas
       *    3. potrf's & trsm's to be batched
       *
       *  The dimensions are stored in lists, in the order specified.
       */

      TIMER_START(tstart,1);      
      /* loop over supernodes */
      for(i = 0; i < nbatch; i++)
      {

        /* get supernode dimensions */
        Int node_ptr = start + node + i;
        Int s 		= supernode_levels[node_ptr];
        Int k1 		= Super [s] ;                          	/* s contains columns k1 to k2-1 of L */
        Int k2 		= Super [s+1] ;
        Int nscol 	= k2 - k1 ;                         	/* # of columns in all of s */
        Int psi 	= Lpi [s] ;                           	/* pointer to first row of s in Ls */
        Int psx 	= Lpx [s] ;                           	/* pointer to first row of s in Lx */
        Int psend 	= Lpi [s+1] ;                   	/* pointer just past last row of s in Ls */
        Int nsrow 	= psend - psi ;     	      		/* # of rows in all of s */
        Int pend 	= psx + nsrow * nscol ;   	        /* s is nsrow-by-nscol */
        Int pk 		= psx ;

        Int m 	= nsrow - nscol ;
        Int n 	= nscol;
        Int lda = nsrow;
        Int ldb = nsrow;

        /* store supernode dimensions */
	super[super_count].s 	 = s;
	super[super_count].k1	 = k1;
        super[super_count].k2 	 = k2;
	super[super_count].psi 	 = psi;
	super[super_count].psx 	 = psx;
	super[super_count].nscol = nscol;
	super[super_count].nsrow = nsrow;
  	super_count++;

        /* store potrf dimensions & pointers */
        potrf[potrf_count].score = n;
        potrf[potrf_count].n 	 = n;
        potrf[potrf_count].lda 	 = lda;
        potrf[potrf_count].A 	 = (double *)(Lx + L_ENTRY*psx);
        potrf_count++;

        /* store trsm dimensions & pointers */
        trsm[trsm_count].score  = n;
        trsm[trsm_count].m 	= m;
        trsm[trsm_count].n 	= n;
        trsm[trsm_count].lda 	= lda;
        trsm[trsm_count].ldb 	= ldb;
        trsm[trsm_count].A 	= (double *)(Lx + L_ENTRY*psx);
        trsm[trsm_count].B 	= (double *)(Lx + L_ENTRY*(psx + n));
        trsm_count++;

	/* increment flops */
        SUM(potrf_flops,1,(double)(n*n*n/3.0));
        SUM(trsm_flops,1,(double)(m*n*n));

      } /* end loop over supernodes */
      TIMER_END(tstart,tend,1);









      /*
       *  Store descendant dimensions:
       *
       *  Stores dimensions of descendants in the following order:
       *    1. syrk's to be sent to cuBlas
       *    2. gemm's to be sent to cuBlas
       *    3. syrk's & gemm's to be batched
       *
       *  First gather all dimensions and sort them by size. Then
       *  loop over ordered list of descendants and store its dimensions
       *  into lists, in the order specified.
       *
       */
#pragma omp critical
{

      TIMER_START(tstart,2);
      /* loop over supernodes */
      for(i = 0; i < nbatch; i++)
      {

        /* get supernode dimensions */
        Int s  		= super[i].s;
        Int k2		= super[i].k2;
        Int nscol       = super[i].nscol;
        Int nsrow2      = super[i].nsrow - super[i].nscol;
        Int psi         = super[i].psi;
        Int dnext 	= Head[s];
        Int ndescendant = ndescendants[s];

        /* loop over descendants */
        for(j=0; j < ndescendant; j++)
        {

          /* get descendant dimensions */
	  Int d 	= dnext;
          Int kd1 	= Super [d] ;                       
          Int kd2 	= Super [d+1] ;
          Int ndcol 	= kd2 - kd1 ;                     
          Int pdi 	= Lpi [d] ;                         
          Int pdx 	= Lpx [d] ;                         
          Int pdend 	= Lpi [d+1] ;                     
          Int ndrow 	= pdend - pdi ;                   

          Int p 	= Lpos [d] ;                          
          Int pdi1 	= pdi + p ;                       
          Int pdx1 	= pdx + p ;                       

          Int pdi2;
          for (pdi2 = pdi1 ; pdi2 < pdend && Ls [pdi2] < k2 ; (pdi2)++) ;
          Int ndrow1	 = pdi2 - pdi1 ;                 
          Int ndrow2 	 = pdend - pdi1 ;                

          /* prepare this supernode d for its next ancestor */
          dnext 	 = Next [d] ;
          Lpos [d] 	 = pdi2 - pdi ;
          if (Lpos [d] < ndrow)
          {
            Int dancestor    = SuperMap [Ls [pdi2]] ;
            Next [d] 	     = Head [dancestor] ;
            Head [dancestor] = d ;
          }

          Int m   = ndrow2-ndrow1;
          Int n   = ndrow1;
          Int k   = ndcol;
          Int lda = ndrow;
          Int ldb = ndrow;
          Int ldc = ndrow2;

          /* store descendant dimensions */
	  desc[desc_count].s 	  = i;
          desc[desc_count].pdi1   = pdi1;
	  desc[desc_count].ndrow1 = ndrow1;
	  desc[desc_count].ndrow2 = ndrow2;
	  desc[desc_count].C 	  = (double *)&h_C[counter];
	  desc_count++;
  
          /* store syrk dimensions & pointers */
	  syrk[syrk_count].score = n*k;
	  syrk[syrk_count].n 	 = n;          
          syrk[syrk_count].k 	 = k;
          syrk[syrk_count].lda   = lda;
          syrk[syrk_count].ldc	 = ldc;          
          syrk[syrk_count].A	 = (double *)(Lx + L_ENTRY*pdx1);
          syrk[syrk_count].C	 = (double *)&h_C[counter];
	  syrk_count++;

          /* store gemm dimensions & pointers */
	  gemm[gemm_count].score = m*n*k;
          gemm[gemm_count].m 	 = m;
          gemm[gemm_count].n 	 = n;
          gemm[gemm_count].k	 = k;
          gemm[gemm_count].lda	 = lda;
          gemm[gemm_count].ldb	 = ldb;
          gemm[gemm_count].ldc	 = ldc;
          gemm[gemm_count].A	 = (double *)(Lx + L_ENTRY*(pdx1 + n));
	  gemm[gemm_count].B	 = (double *)(Lx + L_ENTRY*pdx1);
          gemm[gemm_count].C	 = (double *)(&h_C[counter] + L_ENTRY*n);
	  gemm_count++;

	  /* increment pointer to C buff */
	  counter += n*ldc;

	  /* increment flops */
          SUM(syrk_flops,1,(double)(n*n*k));
          SUM(gemm_flops,1,(double)(2*m*n*k));

        } /* end loop over descendants */


        /* prepare for next supernode */
        if(nsrow2 > 0) {
          Lpos [s]        = nscol ;
          Int sparent     = SuperMap [Ls [psi + nscol]] ;
          Next [s]        = Head [sparent] ;
          Head [sparent]  = s ;
        }

        Head [s] = EMPTY ;



      } /* end loop over supernodes */
      TIMER_END(tstart,tend,2);

} /* end pragma omp critical */

      /* sort lists */
/*    qsort ( syrk, syrk_count, sizeof(struct cholmod_syrk_t), (__compar_fn_t) CHOLMOD(sort_syrk) );*/
/*    qsort ( gemm, gemm_count, sizeof(struct cholmod_gemm_t), (__compar_fn_t) CHOLMOD(sort_gemm) );*/
/*    qsort ( potrf, potrf_count, sizeof(struct cholmod_potrf_t), (__compar_fn_t) CHOLMOD(sort_potrf) );*/
/*    qsort ( trsm, trsm_count, sizeof(struct cholmod_trsm_t), (__compar_fn_t) CHOLMOD(sort_trsm) );*/
    









      /*
       *  Initialize Supernode
       *
       *  Initializes the supernode with the following steps:
       *
       *  1. create Map for each supernode (in the batch)
       *  2. initialize Lx for each supernode
       *
       */

      TIMER_START(tstart,3);
      /* loop over supernodes */
      #pragma omp parallel for num_threads(numThreads) if(nbatch >= numThreads) 
      for(i = 0; i < nbatch; i++)
      {

	Int k, ii, pf, pk, p, q, pend, pfend, imap;

        /* get supernode dimensions */
        Int s 		= super[i].s;
        Int k1 		= super[i].k1;
        Int k2 		= super[i].k2;
        Int psi	 	= super[i].psi;
        Int psx 	= super[i].psx;
        Int nsrow 	= super[i].nsrow;     
        Int *Map 	= &h_Map[i*n];

   
        /* construct Map */
        #pragma omp parallel for num_threads(numThreads) if ( nsrow >= 128 && nbatch < numThreads ) 
        for (k = 0 ; k < nsrow ; k++)
          Map [Ls [psi + k]] = k ;


        /* initialize Lx */
        #pragma omp parallel for private ( p, pend, pfend, pf, ii, j, imap, q ) num_threads(numThreads) if ( k2-k1 >= 64 && nbatch < numThreads ) 
        for (k = k1 ; k < k2 ; k++)
        {
          if (stype != 0)
          {
            p = Ap [k] ;
            pend = (Apacked) ? (Ap [k+1]) : (p + Anz [k]) ;
            for ( ; p < pend ; p++)
            {
              ii = Ai [p] ;
              if (ii >= k)
              {
                imap = Map [ii] ;							
                if (imap >= 0 && imap < nsrow)
                {
                  L_ASSIGN (Lx,(imap+(psx+(k-k1)*nsrow)), Ax,Az,p);
                }
              }
            }
          }
          else
          {
            double fjk[2];
            pf = Fp [k] ;
            pfend = (Fpacked) ? (Fp [k+1]) : (p + Fnz [k]) ;
            for ( ; pf < pfend ; pf++)
            {
              j = Fi [pf] ;
              L_ASSIGN (fjk,0, Fx,Fz,pf) ;					
              p = Ap [j] ;
              pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
              for ( ; p < pend ; p++)
              {
                ii = Ai [p] ;
                if (ii >= k)
                {
                  imap = Map [ii] ;
                  if (imap >= 0 && imap < nsrow)
                  {
                    L_MULTADD (Lx,(imap+(psx+(k-k1)*nsrow)),Ax,Az,p, fjk) ;	
                  }
                }
              }
            }
          }
        }  

  
        /* add beta */
        if (beta [0] != 0.0)
        {
          pk = psx ;
          for (k = k1 ; k < k2 ; k++)
          {
            L_ASSEMBLE (Lx,pk, beta) ;		
            pk += nsrow + 1 ;       		
          }
        }

      } /* end loop over supernodes */
      TIMER_END(tstart,tend,3);









      /*
       *  DSYRK
       *
       *   Perform dsyrk on batch of descendants
       *
       */    
 
      TIMER_START1(tstart1);
      /* omp parallel for small matrices */
#ifdef MKLROOT
      mkl_set_num_threads_local(1);
#else
      openblas_set_num_threads(1);
#endif

      /* loop over syrk's */
      #pragma omp parallel for num_threads(numThreads)  
      for(i = 0; i < syrk_count; i++)
      {

        /* get syrk dimensions */
        Int n 	= syrk[i].n;
        Int k 	= syrk[i].k;
        Int lda = syrk[i].lda;
        Int ldc = syrk[i].ldc;

        double *A = (double *)syrk[i].A;
        double *C = (double *)syrk[i].C;

        double one[2]  = {1.0, 0.0};
	double zero[2] = {0.0, 0.0};

        
        if(n < 256 && k < 256) 
        {       
#ifdef REAL
          BLAS_dsyrk ("L", "N",
          	      n, k,              
        	      one,                        
        	      A, lda,   
        	      zero,                       
        	      C, ldc) ;                
#else
          BLAS_zherk ("L", "N",
          	      n, k,              
        	      one,                        
        	      A, lda,   
        	      zero,                       
        	      C, ldc) ;                
#endif
        }
      } /* end loop over syrk's */




#ifdef MKLROOT
      /* no omp parallel for large matrices */
      mkl_set_num_threads_local(numThreads);
#else
      openblas_set_num_threads(numThreads);
#endif

      /* loop over syrk's */
      for(i = 0; i < syrk_count; i++)
      {

        /* get syrk dimensions */
        Int n 	= syrk[i].n;
        Int k 	= syrk[i].k;
        Int lda = syrk[i].lda;
        Int ldc = syrk[i].ldc;

        double *A = (double *)syrk[i].A;
        double *C = (double *)syrk[i].C;

        double one[2]  = {1.0, 0.0};
	double zero[2] = {0.0, 0.0};


        if(n >= 256 || k >= 256)
        {
#ifdef REAL
          BLAS_dsyrk ("L", "N",
                      n, k,
                      one,
                      A, lda,
                      zero,
                      C, ldc) ;
#else
          BLAS_zherk ("L", "N",
                      n, k,
                      one,
                      A, lda,
                      zero,
                      C, ldc) ;       
#endif
        }
      } /* end loop over syrk's */
      TIMER_END1(tstart1,syrk_time,0);
      TIMER_END1(tstart1,syrk_time,1);










      /*
       *  DGEMM
       *
       *  Perform dgemm on batch of descendants
       *
       */
      
      TIMER_START1(tstart1);

#ifdef MKLROOT
      /* omp parallel for small matrices */
      mkl_set_num_threads_local(1);
#else
      openblas_set_num_threads(1);
#endif

      /* loop over gemm's */
      #pragma omp parallel for num_threads(numThreads) 
      for(i = 0; i < gemm_count; i++)
      {

        /* get gemm dimensions */
        Int m 	= gemm[i].m;
        Int n 	= gemm[i].n;
        Int k 	= gemm[i].k;
        Int lda = gemm[i].lda;
        Int ldb = gemm[i].ldb;
        Int ldc = gemm[i].ldc;

        double *A = (double *)gemm[i].A;
        double *B = (double *)gemm[i].B;
        double *C = (double *)gemm[i].C;

        double one[2] = {1.0, 0.0};
	double zero[2] = {0.0, 0.0};


	if(m*n*k < 256*256*256) 
	{
          if (m > 0) 
          {
#ifdef REAL
            BLAS_dgemm ("N", "C",
          	        m, n, k,          
          	        one,                            
          	        A, lda,                          
          	        B, ldb,                          
          	        zero,                           
          	        C, ldc) ;
#else
            BLAS_zgemm ("N", "C",
                        m, n, k,          
                        one,                            
                        A, lda,                          
                        B, ldb,
                        zero,                           
                        C, ldc) ;
#endif
          }
        }
      } /* end loop over gemm's */




#ifdef MKLROOT
      /* no omp parallel for large matrices */
      mkl_set_num_threads_local(numThreads);
#else
      openblas_set_num_threads(numThreads);
#endif

      /* loop over gemm's */
      for(i = 0; i < gemm_count; i++)
      {

        /* get gemm dimensions */
        Int m 	= gemm[i].m;
        Int n 	= gemm[i].n;
        Int k 	= gemm[i].k;
        Int lda = gemm[i].lda;
        Int ldb = gemm[i].ldb;
        Int ldc = gemm[i].ldc;

        double *A = (double *)gemm[i].A;
        double *B = (double *)gemm[i].B;
        double *C = (double *)gemm[i].C;

        double one[2] = {1.0, 0.0};
	double zero[2] = {0.0, 0.0};


	if(m*n*k >= 256*256*256) 
	{
          if (m > 0)
          {
#ifdef REAL
            BLAS_dgemm ("N", "T",
                        m, n, k,
                        one,
                        A, lda,
                        B, ldb,
                        zero,
                        C, ldc) ;
#else
            BLAS_zgemm ("N", "C",
                        m, n, k,
                        one,
                        A, lda,
                        B, ldb,
                        zero,
                        C, ldc) ;
#endif
          }
        }
      } /* end loop over gemm */
      TIMER_END1(tstart1,gemm_time,0);
      TIMER_END1(tstart1,gemm_time,1);









      /*
       *  Assembly
       *
       *  Assemble schur complements of a batch of descendants
       *
       */

      TIMER_START(tstart,4);
       /* loop over descendants */
      for(i = 0; i < desc_count; i++)      
      {

        Int ii, j, q, px, ix;

        /* get supernode dimensions */
        Int s = desc[i].s;
        Int psx = super[s].psx;
        Int nsrow = super[s].nsrow;
        Int *Map = &h_Map[s*n];

        /* get descendant dimensions */
        Int pdi1 = desc[i].pdi1;
        Int ndrow1 = desc[i].ndrow1;
        Int ndrow2 = desc[i].ndrow2;

        double *C = (double *)desc[i].C;


        #pragma omp parallel for private ( j, ii, px, q ) num_threads(numThreads) if (ndrow1 > 64 )
        for (j = 0 ; j < ndrow1 ; j++)                      
        {
          px = psx + Map [Ls [pdi1 + j]]*nsrow ;
          for (ii = j ; ii < ndrow2 ; ii++)                   
          {
            q = px +  Map [Ls [pdi1 + ii]] ;
            ix = ii+ndrow2*j;
            Lx[q] -= C[ix];
          }
        }

      } /* end loop over descendants */
      TIMER_END(tstart,tend,4);









      /*
       *  Cholesky Factorization
       *
       *  Perform dpotrf for a batch of supernodes.
       *
       */

      TIMER_START1(tstart1);

#ifdef MKLROOT
      /* omp parallel for small matrices */
      mkl_set_num_threads_local(1);
#else
      openblas_set_num_threads(1);
#endif

      /* loop over potrf's */
      #pragma omp parallel for num_threads(numThreads)
      for(i = 0; i < potrf_count; i++)
      {

        /* get potrf dimensions */        
        Int n 	  = potrf[i].n;
        Int lda   = potrf[i].lda;      

        double *A = (double *)potrf[i].A;   
       
 
        if(potrf_count >= numThreads) 
	{
#ifdef REAL
          LAPACK_dpotrf ("L",
                         n,
                         A, lda,
                         iinfo[i]) ;

#else
          LAPACK_dpotrf ("L",
                         n,
                         A, lda,
                         iinfo[i]) ;
#endif
        }
      } /* end loop over potrf's */



#ifdef MKLROOT
      /* no omp parallel for large matrices */
      mkl_set_num_threads_local(numThreads);
#else
      openblas_set_num_threads(numThreads);
#endif

      /* loop over potrf's */
      for(i = 0; i < potrf_count; i++)
      {

        /* get potrf dimensions */
        Int n 	  = potrf[i].n;
        Int lda   = potrf[i].lda;

        double *A = (double *)potrf[i].A;


        if(potrf_count < numThreads) 
	{
#ifdef REAL
          LAPACK_dpotrf ("L",
                         n,
                         A, lda,
                         iinfo[i]) ;

#else
          LAPACK_dpotrf ("L",
                         n,
                         A, lda,
                         iinfo[i]) ;
#endif
        }
      } /* end loop potrf's */
      TIMER_END1(tstart1,potrf_time,0);
      TIMER_END1(tstart1,potrf_time,1);    









      /*
       *  Triangular Solve
       *
       *  Perform dtrsm for a batch of supernodes.
       *
       */

      TIMER_START1(tstart1);

#ifdef MKLROOT
      /* omp parallel for small matrices */
      mkl_set_num_threads_local(1);
#else
      openblas_set_num_threads(1);
#endif

      /* loop over trsm's */
      #pragma omp parallel for num_threads(numThreads)  
      for(i = 0; i < trsm_count; i++)
      {

        /* get trsm dimensions */
        Int m 	= trsm[i].m;
        Int n	= trsm[i].n;
        Int lda = trsm[i].lda;
        Int ldb = trsm[i].ldb;

        double *A = (double *)trsm[i].A;
        double *B = (double *)trsm[i].B;

        double one[2] = {1.0, 0.0};


        if(trsm_count >= numThreads && m*n < 256*256)         
        {
          if (m > 0) 
          {
#ifdef REAL
             BLAS_dtrsm ("R", "L", "T", "N",
                         m, n,                   
                         one,                               
                         A, lda,         
                         B, ldb) ;
#else
             BLAS_ztrsm ("R", "L", "C", "N",
                         m, n,                   
                         one,                               
                         A, lda,         
                         B, ldb) ;
#endif
          }
        }
      } /* end loop over trsm's */



#ifdef MKLROOT
      /* no omp parallel for large matrices */
      mkl_set_num_threads_local(numThreads);
#else
      openblas_set_num_threads(numThreads);
#endif

      /* loop over trsm's */
      for(i = 0; i < trsm_count; i++)
      {

        /* get trsm dimensions */
        Int m 	= trsm[i].m;
        Int n 	= trsm[i].n;
        Int lda = trsm[i].lda;
        Int ldb = trsm[i].ldb;

        double *A = (double *)trsm[i].A;
        double *B = (double *)trsm[i].B;

        double one[2] = {1.0, 0.0};


        if(trsm_count < numThreads || m*n >= 256*256) 
	{
          if (m > 0)
          {
#ifdef REAL
             BLAS_dtrsm ("R", "L", "T", "N",
                         m, n,
                         one,
                         A, lda,
                         B, ldb) ;
#else
             BLAS_ztrsm ("R", "L", "C", "N",
                         m, n,
                         one,
                         A, lda,
                         B, ldb) ;
#endif
          }
        }
      } /* end loop over trsm's */
      TIMER_END1(tstart1,trsm_time,0);
      TIMER_END1(tstart1,trsm_time,1);










      /*
       *  Prepare for next batch of supernodes
       */

      /* loop over supernodes */
      for(i = 0; i < nbatch; i++)
      {

        /* get supernode dimensions */
        Int s     	= super[i].s;
        Int k1		= super[i].k1;
        Int nscol 	= super[i].nscol;
        Int nsrow2      = super[i].nsrow - super[i].nscol;
        Int psi   	= super[i].psi;
/*
        if(nsrow2 > 0) {
          Lpos [s]	  = nscol ;
          Int sparent 	  = SuperMap [Ls [psi + nscol]] ;
          Next [s] 	  = Head [sparent] ;
          Head [sparent]  = s ;
        }

        Head [s] = EMPTY ;  
*/        

        /* check if the matrix is not positive definite */
        if (iinfo[i] != 0)
        {
          /* Matrix is not positive definite.  dpotrf/zpotrf do NOT report an
           * error if the diagonal of L has NaN's, only if it has a zero. */
          L->minor = k1 + iinfo[i] - 1;
          ERROR (CHOLMOD_NOT_POSDEF, "matrix not positive definite") ;
          return Common->status;/*(Common->status >= CHOLMOD_OK)*/
        }

      } /* end loop over supernodes */










      /* increment supernode id */
      node += nbatch;
      TIMER_END(tstart,tend,0);      
 
    } /* end loop over batches*/




    /* increment flops */
    SUM(syrk_flops,0,(double)(syrk_flops[1]));
    SUM(gemm_flops,0,(double)(gemm_flops[1]));
    SUM(potrf_flops,0,(double)(potrf_flops[1]));
    SUM(trsm_flops,0,(double)(trsm_flops[1]));

    /* print benchmarks per level */
    PRINTFV("%d\t",level);    
    PRINTFV("%d\t",nbatch);
    PRINTFV("%d\t",end-start);
    PRINTFV("%f\t",syrk_time[1]);
    PRINTFV("%f\t",gemm_time[1]);
    PRINTFV("%f\t",potrf_time[1]);
    PRINTFV("%f\t",trsm_time[1]);
    PRINTFV("%f\t",1.0e-9*syrk_flops[1]/syrk_time[1]);
    PRINTFV("%f\t",1.0e-9*gemm_flops[1]/gemm_time[1]);
    PRINTFV("%f\t",1.0e-9*potrf_flops[1]/potrf_time[1]);
    PRINTFV("%f\n",1.0e-9*trsm_flops[1]/trsm_time[1]);

    



    /* clear timers */
    CLEAR1(syrk_time,1);
    CLEAR1(gemm_time,1);
    CLEAR1(potrf_time,1);
    CLEAR1(trsm_time,1);

  } /* end loop over levels */




  /* print overall benchmarks */
  PRINTFV("\n\nSubtree %d benchmarks:\n",subtree);
  PRINTF("\n- time -\n");
  PRINTFV("total:        %f\n",tend[0]);
  PRINTFV("superDim:     %f\n",tend[1]);
  PRINTFV("descDim:      %f\n",tend[2]);
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
  return Common->status; /*(Common->status >= CHOLMOD_OK) ;*/
}


/*
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
*/











