/* ========================================================================== */
/* === Supernodal/t_cholmod_super_numeric =================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Supernodal Module.  Copyright (C) 2005-2012, Timothy A. Davis
 * The CHOLMOD/Supernodal Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * ---------------------------------------------------------------------------*/


/*
 *
 * Description:
 *   Contains functions for factorization
 *   of the elimination tree
 *
 */


/* include */

#include "cholmod_template.h"











/*
 * Function:
 *   cholmod_super_numeric
 *
 * Description:
 *   Factorizes elimination tree in one of two ways:
 *   1. Splits tree into subtree and:
 *      a. factorize subtree with GPU subtree algorithm
 *      b. factorize top-of-tree subtree with root algorithm
 *   2. Factorizes entire tree with CPU algorithm. 
 *
 */
static int TEMPLATE (cholmod_super_numeric)
(
 cholmod_sparse *A,  			/* matrix to factorize */
 cholmod_sparse *F,  			/* F = A' or A(:,f)' */
 double beta [2],    			/* beta*I is added to diagonal of matrix to factorize */
 cholmod_factor *L,  			/* factorization */
 cholmod_dense *Cwork,       		/* size (L->maxcsize)-by-1 */
 cholmod_common *Common
 )
{





  /* global variables */
  Int i, j, k, size ; 
  Int *LpxSub, *Iwork;
  struct cholmod_subtree_order_t *Bwork;
  double *tstart, *tend, *bstart, *bend, *Xwork;
  struct cholmod_global_pointers *gb_p, gb_pointer_struct;
  struct cholmod_cpu_pointers *cpu_p, cpu_pointer_struct;
  struct cholmod_gpu_pointers *gpu_p, gpu_pointer_struct;
  struct cholmod_tree_pointers *tree_p, tree_pointer_struct;
  struct cholmod_profile_pointers *prof_p, prof_pointer_struct;
  struct cholmod_loadbalance_pointers *lb_p, lb_pointer_struct;





  /* set structure pointers */
  gb_p   	= &gb_pointer_struct ;
  cpu_p  	= &cpu_pointer_struct ;
  gpu_p  	= &gpu_pointer_struct ;
  tree_p 	= &tree_pointer_struct ;
  prof_p 	= &prof_pointer_struct ;
  lb_p   	= &lb_pointer_struct ;

  /* clear global variables */
  gb_p->runType   = 0;  	
  gb_p->numGPU    = 0;
  gb_p->numDevice = 0;
  gb_p->numSubtree = 0;
  gb_p->numRoot   = 0;
  gb_p->work_size = 0;

  gb_p->maxCsize  = 0;
  gb_p->maxndesc  = 0;
  gb_p->maxbatch  = 0;
  gb_p->maxnsrow  = 0;
  gb_p->maxnscol  = 0;

  for(i=0; i < CHOLMOD_MAX_NUM_GPUS; i++) gb_p->check[i] = 0;




  PRINTF("\n\n\n");
  PRINTFV("useGPU: %d\n",Common->useGPU);
  PRINTFV("numGPU: %d\n",Common->numGPU);
  PRINTFV("useHybrid: %d\n",Common->useHybrid);
  PRINTFV("ompNumThreads: %d\n",Common->ompNumThreads);
  PRINTFV("partialFactorization: %d\n",Common->partialFactorization);
  PRINTFV("maxGpuMemBytes: %ld\n",Common->maxGpuMemBytes);
  
  /* hybrid is enabled */
  if(Common->useHybrid == 1) {  
    gb_p->runType = 0;          /* set to hybrid */
  }
  else {
    gb_p->runType = 2;          /* set to GPU only */    
  }

  /* not enough supernodes in the elimination tree */
  if(L->nsuper <= SUPERNODE_MIN) {
    gb_p->runType = -1;          /* set to CPU serial */    
  }

  /* GPU is not enabled */
  if(Common->numGPU == 0 || Common->useGPU == 0) {
    gb_p->runType = 1;          /* set to CPU only */    
  }

  /* matrix is complex */
  #ifdef COMPLEX
    if(gb_p->runType != 1 && gb_p->runType != -1) {    
      gb_p->runType = 3;        /* set to root only */      
    }
  #endif

  /* determine whether to use CPU serial */
  if((Common->ompNumThreads == 1 && Common->useGPU == 0) || Common->partialFactorization == 1) {
    gb_p->runType = -1;		/* set to CPU serial */    
  }

  /* GPU is not enabled */
  #ifndef SUITESPARSE_CUDA
    if(Common->partialFactorization == 1)
      gb_p->runType = -1;
    else
      gb_p->runType = 1;
  #endif


  


  /* print type of run */
  PRINTFV("\nrunType:%d\t",gb_p->runType);
  if(gb_p->runType == 0)      PRINTF("GPU + CPU (hybrid)\n");
  if(gb_p->runType == 1)      PRINTF("CPU only\n");
  if(gb_p->runType == 2)      PRINTF("GPU only\n");
  if(gb_p->runType == 3)      PRINTF("root only\n");





  /* allocate memory for subtree algorithm */
  if(gb_p->runType != -1) {  		/* only if subtree algorithm chosen */

    /* determine size for load-balance arrays*/
    if(L->nsuper < CHOLMOD_MAX_NUM_GPUS) size = CHOLMOD_MAX_NUM_GPUS+2;
    else				 size = L->nsuper;

    /* allocate workspace */
    gb_p->IworkSize                       = 25*(L->nsuper + 1) + (Common->numGPU+4)*(size + 1);
    gb_p->XworkSize                       = 2*(L->nsuper + 1) + (size + 1);
    gb_p->BworkSize                       = L->nsuper;

    gb_p->Iwork                           = CHOLMOD(malloc) (gb_p->IworkSize, sizeof (Int), Common) ;
    gb_p->Xwork                           = CHOLMOD(malloc) (gb_p->XworkSize, sizeof (double), Common) ;
    gb_p->Bwork                           = CHOLMOD(malloc) (gb_p->BworkSize, sizeof (struct cholmod_subtree_order_t), Common) ;    

    Iwork                                 = gb_p->Iwork;
    Xwork				  = gb_p->Xwork;
    Bwork				  = gb_p->Bwork;

    /* check if enough memory */
    if (Common->status < CHOLMOD_OK)
    {      
      gb_p->Iwork = CHOLMOD(free) (gb_p->IworkSize, sizeof (Int), gb_p->Iwork, Common) ;
      gb_p->Xwork = CHOLMOD(free) (gb_p->XworkSize, sizeof (double), gb_p->Xwork, Common) ;
      gb_p->Bwork = CHOLMOD(free) (gb_p->BworkSize, sizeof (struct cholmod_subtree_order_t), gb_p->Bwork, Common) ;
      return (FALSE) ;
    }


    /* clear workspace */
    memset(Iwork,0,gb_p->IworkSize*sizeof(Int));
    memset(Xwork,0,gb_p->XworkSize*sizeof(double));
    memset(Bwork,0,gb_p->BworkSize*sizeof(struct cholmod_subtree_order_t));


    tree_p->supernode_subtree             = Iwork;
    tree_p->supernode_subtree_ptrs        = Iwork + 1*(size_t)(L->nsuper + 1);
    tree_p->supernode_batch               = Iwork + 2*(size_t)(L->nsuper + 1);
    tree_p->supernode_levels              = Iwork + 3*(size_t)(L->nsuper + 1);
    tree_p->supernode_levels_ptrs         = Iwork + 4*(size_t)(L->nsuper + 1);
    tree_p->supernode_levels_subtree_ptrs = Iwork + 5*(size_t)(L->nsuper + 1);
    tree_p->supernode_parent              = Iwork + 6*(size_t)(L->nsuper + 1);
    tree_p->supernode_children            = Iwork + 7*(size_t)(L->nsuper + 1);
    tree_p->supernode_children_ptrs       = Iwork + 8*(size_t)(L->nsuper + 1);
    tree_p->supernode_children_num        = Iwork + 9*(size_t)(L->nsuper + 1);
    tree_p->supernode_children_num2       = Iwork + 10*(size_t)(L->nsuper + 1);
    tree_p->supernode_children_count      = Iwork + 11*(size_t)(L->nsuper + 1);
    tree_p->supernode_children_count2     = Iwork + 12*(size_t)(L->nsuper + 1);
    tree_p->supernode_num_levels          = Iwork + 13*(size_t)(L->nsuper + 1);
    tree_p->level_descendants             = Iwork + 14*(size_t)(L->nsuper + 1);
    tree_p->level_descendants_ptrs        = Iwork + 15*(size_t)(L->nsuper + 1);
    tree_p->level_num_desc                = Iwork + 16*(size_t)(L->nsuper + 1);
    tree_p->level_num_desc_ptrs           = Iwork + 17*(size_t)(L->nsuper + 1);
    tree_p->supernode_size_desc           = Iwork + 18*(size_t)(L->nsuper + 1);
    tree_p->supernode_size                = Iwork + 19*(size_t)(L->nsuper + 1);
    tree_p->supernode_root                = Iwork + 20*(size_t)(L->nsuper + 1);
    tree_p->factor_size                   = Iwork + 21*(size_t)(L->nsuper + 1);
    tree_p->ndescendants                  = Iwork + 22*(size_t)(L->nsuper + 1);
    lb_p->numSubtreePerDevice             = Iwork + 23*(size_t)(L->nsuper + 1);
    lb_p->listSubtreePerDevice            = Iwork + 23*(size_t)(L->nsuper + 1) + (size_t)(size + 1);
    LpxSub                                = Iwork + 23*(size_t)(L->nsuper + 1) + (Common->numGPU+4)*(size_t)(size + 1);;

    tree_p->supernode_flop                = Xwork;
    lb_p->subtreeSize                     = Xwork + (size_t)(L->nsuper + 1);
    lb_p->workPerDevice                   = Xwork + 2*(size_t)(L->nsuper + 1);

    lb_p->subtreeReorder                  = Bwork;
    
  }





  /* allocate integer workspace */
  Iwork        	      = Common->Iwork;
  cpu_p->SuperMap     = Iwork;                                     
  cpu_p->RelativeMap  = Iwork + L->n;
  cpu_p->Next         = Iwork + 2*((size_t)L->n);                     
  cpu_p->Lpos         = Iwork + 2*((size_t)L->n) + L->nsuper;            
  cpu_p->Next_save    = Iwork + 2*((size_t)L->n) + 2*((size_t)L->nsuper);
  cpu_p->Lpos_save    = Iwork + 2*((size_t)L->n) + 3*((size_t)L->nsuper);
  cpu_p->Previous     = Iwork + 2*((size_t)L->n) + 4*((size_t)L->nsuper);





  /* set host pointers */    
  cpu_p->C 	= Cwork->x ;
  cpu_p->Map  	= Common->Flag ;   
  cpu_p->Head 	= Common->Head ;   
  cpu_p->Ls 	= L->s ;
  cpu_p->Lpi 	= L->pi ;
  cpu_p->Lpx 	= LpxSub;
  cpu_p->Super 	= L->super ;
  cpu_p->Lx 	= L->x ;
  cpu_p->stype 	= A->stype ;
  cpu_p->beta 	= beta;

  cpu_p->Ap             = A->p ;
  cpu_p->Ai             = A->i ;
  cpu_p->Ax             = A->x ;
  cpu_p->Az             = A->z ;
  cpu_p->Anz            = A->nz ;
  cpu_p->Apacked        = A->packed ;

  if (cpu_p->stype != 0)
  {
      cpu_p->Fp 	= NULL ;
      cpu_p->Fi 	= NULL ;
      cpu_p->Fx 	= NULL ;
      cpu_p->Fz 	= NULL ;
      cpu_p->Fnz 	= NULL ;
      cpu_p->Fpacked 	= TRUE ;
  }
  else
  {
      cpu_p->Fp 	= F->p ;
      cpu_p->Fi 	= F->i ;
      cpu_p->Fx 	= F->x ;
      cpu_p->Fz 	= F->z ;
      cpu_p->Fnz 	= F->nz ;
      cpu_p->Fpacked 	= F->packed ;
  }





  /* set timer pointers */
  tstart        = prof_p->g_start;
  tend          = prof_p->g_end;
  bstart        = prof_p->b_start;
  bend          = prof_p->b_end;





  /* check if functionality available - (not supported for GPU subtree) */
  if(cpu_p->Apacked==0 || cpu_p->stype==0 || cpu_p->beta[0]!=0) {
    if(gb_p->runType != 1 && gb_p->runType != -1) {    
      gb_p->runType = 3;               			 /* set to root only */      
    }   
  }





#ifdef SUITESPARSE_CUDA
  /* clear floating point exceptions */
  if (feclearexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_DIVBYZERO | FE_INVALID | FE_INEXACT | FE_ALL_EXCEPT)){
    PRINTF("\nfloating-point exceptions not cleared!\n");
  }
  else{
    PRINTF("\nfloating-point exceptions cleared!\n");
  }
#endif





  /* clear the Map so that changes in the pattern of A can be detected */
  #pragma omp parallel for num_threads(Common->ompNumThreads) if ( L->n > 128 ) schedule (static)
  for (i = 0 ; i < L->n ; i++)
    cpu_p->Map [i] = EMPTY ;










  /*
   * Serial Factorization
   *
   * Description:
   * Performs serial factorization on the elimination tree.
   * Steps:
   *   1. factorize elimination tree serially
   */
  if(gb_p->runType == -1)  
  {
    PRINTF("\n\n\nSERIAL FACTORIZATION selected..\n");
    int deviceid = 0, check = 0;
    check = TEMPLATE2 (CHOLMOD(gpu_factorize_cpu_serial))( Common, L, gb_p, cpu_p, tree_p, prof_p, deviceid);
    if(check) return (Common->status >= CHOLMOD_OK);	/* early exit if not positive-definite */
  }










  /*
   * Parallel Factorization
   *
   * Description:
   * Performs parallel factorization on the elimination tree.
   * Steps:
   *   1. build elimination tree
   *   2. build subtrees (through binary search)
   *   3. load balance devices
   *   4. initialize CPU & GPU
   *   5. factorize subtrees in parallel
   *   6. factorize root
   */
  if(gb_p->runType != -1)  
  {

    PRINTF("\n\n\nPARALLEL FACTORIZATION selected..\n");
    /* start factorize timer.. */
    TIMER_START(tstart,0);


    /* 
     * Build elimination tree
     *
     * Description:
     *   stores information about elimination tree:
     *   supernode sizes, # descendants, # children, children, parents, root.
     */
    PRINTF("\n\n\nbuild elimination tree..\n");
    TIMER_START(tstart,1);
    TEMPLATE2 (CHOLMOD (build_tree))( Common,L,gb_p,cpu_p,tree_p );

    /* store copy of # children per supernode */   
    memcpy(tree_p->supernode_children_num2, tree_p->supernode_children_num, L->nsuper*sizeof(Int));
    TIMER_END(tstart,tend,1);










    /* 
     * Binary search for optimal subtree size
     *  
     * Description:
     * perform binary search to find optimal subtree size. Performs up to BINARY_SEARCH 
     * steps.
     * 
     */ 
    PRINTF("\n\n\nprocess subtree (binary search) ..\n"); 
    TIMER_START(tstart,2);  
    TEMPLATE2 (CHOLMOD(binarysearch_tree))( Common, A, L, gb_p, cpu_p, tree_p, LpxSub);
    TIMER_END(tstart,tend,2);  










    /* 
     * Load-balance Devices
     *
     * Description:
     * Reorder subtree (subtrees) by size, which is quantified by its workload (flop/flops).
     * Then load-balance subtree to different device (GPU & CPU), for maximum utilization.   
     */
    PRINTF("\n\n\nload-balance devices..\n");

    TIMER_START(tstart,3);  
    TEMPLATE2 (CHOLMOD(loadbalance_gpu))( Common,gb_p,tree_p,lb_p);
    TIMER_END(tstart,tend,3);  










    /*
     * Initialize GPU & CPU
     *
     * Description:
     * 1. initialize GPU (set pointers, copy memory, etc.)
     * 2. initialize CPU (clear Lx factor, allocate memory for parallel CPU algorithm)
     */    
    PRINTF("\n\n\ninit GPU & CPU..\n");
    TIMER_START(tstart,4);  
    TEMPLATE2 (CHOLMOD(initialize_gpu))(Common,L,A,gb_p,gpu_p,cpu_p);	/* initialize GPU */
    TEMPLATE2 (CHOLMOD(initialize_cpu))(Common,L,gb_p,cpu_p,tree_p);	/* initialize CPU */
    TIMER_END(tstart,tend,4);  










    /* print system information */
    PRINTF("\n\n\nfactorize tree..\n");
    PRINTFV("total # supernodes: %d\n",L->nsuper);
    PRINTFV("numSubtree: %d\n",gb_p->numSubtree);
    PRINTFV("numDevice:	%d\n",gb_p->numDevice);
    for(i = 0; i < Common->numGPU+2; i++) {
      PRINTFV("device:%d ",i);
      PRINTFV("numSubtreePerDevice:%d ",lb_p->numSubtreePerDevice[i]);
      PRINTFV("workPerDevice:%d\n",lb_p->workPerDevice[i]);
    }

    PRINTF("\n\ntype of run: ");
    if(gb_p->runType == 0)	PRINTF("GPU + CPU (hybrid)\n");
    if(gb_p->runType == 1)      PRINTF("CPU only\n");    
    if(gb_p->runType == 2)      PRINTF("GPU only\n");
    if(gb_p->runType == 3)      PRINTF("root only\n");    











    /* 
     * Supernodal numerical factorization (with GPU & CPU)
     *  
     * Description:
     * factorization using three algorithms:
     * 1. GPU (subtree that fits GPU)
     * 2. CPU (subtree concurrent with GPU)
     * 3. root (CPU/GPU) (last subtree that does not fit GPU)       
     *
     * If root_only or CPU_only = 1, the factorization is done
     * entirely on the root or CPU.  
     */
    /* start timer for factorization */
    PRINTF("\n\n\nsupernodal numerical factorization..\n");
    TIMER_START(tstart,5);  

    omp_set_nested(1);     		/* allow for nested omp */

    /* set # omp threads: 
     *   1. CPU only:     1
     *   2. GPU only:     Common->numGPU 
     *   3. hybrid:       Common->numGPU + 1  
     */
    if(gb_p->runType == 1)      gb_p->numDevice = 1;                            /* CPU only */ 
    else if(gb_p->runType == 2) gb_p->numDevice = Common->numGPU;               /* GPU only */ 
    else              		gb_p->numDevice = Common->numGPU + 1;		/* GPU + CPU (hybrid) */





    /* loop over all devices (GPU,CPU) */
    #pragma omp parallel num_threads(gb_p->numDevice)
    {        
      /* local variables */
      int deviceid, subtreeid, numSubtreePerDevice, check = 0;

      /* set variables */
      deviceid = omp_get_thread_num();				/* set device id*/
      gpu_p->gpuid = deviceid;					/* gpuid (for GPU algorithm) */
      numSubtreePerDevice = (int)(lb_p->numSubtreePerDevice[deviceid]);




      /*
       * GPU subtree algorithm
       *
       * Description:
       * Performs factorization on subtree of the elimination tree.
       * Uses GPU only algorithm. Optimized for small matrices.
       * Case where subtree of elimination tree fits the GPU. Is
       * optimized for small matrices.
       *
       */
      if(deviceid < Common->numGPU)      
      {

        /* set device */
#ifdef SUITESPARSE_CUDA
        cudaSetDevice(deviceid);
#endif

        /* loop over subtree in current GPU device */
        for(subtreeid = 0; subtreeid < numSubtreePerDevice; subtreeid++)        
        {
          /* get current subtree & # supernodes */
          Int subtree 	= lb_p->listSubtreePerDevice[subtreeid + deviceid*gb_p->numSubtree];
          Int numSuper 	= tree_p->supernode_subtree_ptrs[subtree+1] - tree_p->supernode_subtree_ptrs[subtree];
  
          PRINTF("\n\nGPU start -\t");
          PRINTFV("device:%d ",deviceid);
          PRINTFV("subtree:%d ",subtree);

          TIMER_START(bstart,deviceid);
          TEMPLATE2 (CHOLMOD(gpu_factorize_subtree))( Common, gb_p, gpu_p, cpu_p, tree_p, prof_p, L, deviceid, numSuper, subtree, LpxSub);          
    	  TIMER_END(bstart,bend,deviceid);        

          PRINTF("\n\nGPU end -\t");
          PRINTFV("device:%d ",deviceid);
          PRINTFV("subtree:%d ",subtree);
          PRINTFV("nsuper:%d ",numSuper);
          PRINTFV("subtreeSize:%f ",lb_p->subtreeSize[subtree]);
          PRINTFV("time:%f\n",bend[deviceid]);
        } /* end loop over subtree */
      } /* end if GPU subtree */





      /*
       * CPU algorithm
       *
       * Description:
       * Performs factorization on subtree of the elimination tree.
       * Uses CPU only algorithm. Goal of utilizing CPU while GPU
       * is busy. If CPU_only = 1, performs factorization on entire
       * tree. 
       *
       * Call one of two functions:
       *   1. gpu_factorize_cpu_serial (serial factorization)
       *   2. gpu_factorize_cpu_parallel (parallel factorization) 
       * 
       */
      if(deviceid == Common->numGPU)      
      {

        /* loop over subtree in CPU device */
        for(subtreeid = 0; subtreeid < numSubtreePerDevice; subtreeid++)
        {
        
          /* get current subtree & # supernodes */
          Int subtree 	= lb_p->listSubtreePerDevice[subtreeid + deviceid*gb_p->numSubtree];
          Int numSuper 	= tree_p->supernode_subtree_ptrs[subtree+1] - tree_p->supernode_subtree_ptrs[subtree];

          PRINTF("\n\nCPU start -\t");
          PRINTFV("device:%d ",deviceid);
          PRINTFV("subtree:%d ",subtree);

    	  TIMER_START(bstart,deviceid);        
          check = TEMPLATE2 (CHOLMOD(gpu_factorize_cpu_parallel))( Common, L, gb_p, cpu_p, tree_p, prof_p, deviceid, subtree);
  	  TIMER_END(bstart,bend,deviceid);        

          PRINTF("\n\nCPU end -\t");
          PRINTFV("device:%d ",deviceid);
          PRINTFV("subtree:%d ",subtree);
          PRINTFV("nsuper:%d ",numSuper);
          PRINTFV("subtreeSize:%f ",lb_p->subtreeSize[subtree]);
          PRINTFV("time:%f\n",bend[deviceid]);

          if(check) gb_p->check[deviceid] = check;
        } /* end loop over subtree */
      } /* end if CPU subtree */
    } /* end loop over devices (OMP threads) */




    /* early exit if subtree not positive-definite */
    for(i=0; i < CHOLMOD_MAX_NUM_GPUS; i++) {
      if(gb_p->check[i]) return (Common->status >= CHOLMOD_OK);
    }





    /*
     * root algorithm
     *
     * Description:
     * Performs factorization on top-of-tree subtree of the
     * elimination tree. Uses CPU/GPU algorithm. Optimized
     * for large matrices. Case where subtree does not fit
     * the GPU. If root_only = 1, performs factorization on
     * entire tree.
     *
     */
    int deviceid = Common->numGPU+1;
    int subtreeid, check = 0;
    int numSubtreePerDevice = (int)(lb_p->numSubtreePerDevice[deviceid]);

    /* reset Cbuff for root algorithm */
    /*cpu_p->C      = Cwork->x ;*/

    if(deviceid == Common->numGPU+1)
    {

      /* wait until all subtree are factorized */
 
      /* loop over subtree in root */
      for(subtreeid = 0; subtreeid < numSubtreePerDevice; subtreeid++)
      {

        /* get current subtree & # supernodes */
        Int subtree 	= lb_p->listSubtreePerDevice[subtreeid + deviceid*gb_p->numSubtree];
        Int numSuper 	= tree_p->supernode_subtree_ptrs[subtree+1] - tree_p->supernode_subtree_ptrs[subtree];
 
        PRINTF("\n\nroot start -\t");
        PRINTFV("device:%d ",deviceid);
        PRINTFV("subtree:%d ",subtree);

        TIMER_START(bstart,deviceid);      
        check = TEMPLATE2 (CHOLMOD(gpu_factorize_root_parallel))( Common, L, gpu_p, cpu_p, tree_p, subtree );
        TIMER_END(bstart,bend,deviceid);      

        PRINTF("\n\nroot end -\t");
        PRINTFV("device:%d ",deviceid);
        PRINTFV("subtree:%d ",subtree);
        PRINTFV("nsuper:%d ",numSuper);
        PRINTFV("subtreeSize:%f ",lb_p->subtreeSize[subtree]);
        PRINTFV("time:%f\n",bend[deviceid]);

        if(check) return (Common->status >= CHOLMOD_OK);     /* early exit if not positive-definite */

      } /* end loop over subtree */
    } /* end if root subtree */





    TIMER_END(tstart,tend,5);
    TIMER_END(tstart,tend,0);





    /* Print timers */
    PRINTF("\n\n\n");
    PRINTFV("total:               \t%f\n",tend[0]);
    PRINTFV("construct tree:      \t%f\n",tend[1]);
    PRINTFV("construct subtree:   \t%f\n",tend[2]);
    PRINTFV("load-balance:        \t%f\n",tend[3]);
    PRINTFV("init GPU & CPU:      \t%f\n",tend[4]);
    PRINTFV("factorize:           \t%f\n",tend[5]);
    PRINTF("\n");

  } /* end if parallel factorization */










  /* success; matrix is positive definite */
  L->minor = L->n ;


  PRINTF("\n\n\nfree GPU & CPU..\n");

#ifdef SUITESPARSE_CUDA
  /* finalize gpu */    
  CHOLMOD (gpu_end) (Common) ;
#endif

  /* free arrays used for subtree algorithm */
  if(gb_p->runType != -1)  
  {
    gb_p->Iwork = CHOLMOD(free) (gb_p->IworkSize, sizeof (Int), gb_p->Iwork, Common) ;
    gb_p->Xwork = CHOLMOD(free) (gb_p->XworkSize, sizeof (double), gb_p->Xwork, Common) ;
    gb_p->Bwork = CHOLMOD(free) (gb_p->BworkSize, sizeof (struct cholmod_subtree_order_t), gb_p->Bwork, Common) ;

/*    if(gb_p->runType != 3 && gb_p->runType != 2) */
    {   
      gb_p->Cwork = CHOLMOD(free) (gb_p->CworkSize, sizeof (double), gb_p->Iwork, Common) ;
      gb_p->Mapwork = CHOLMOD(free) (gb_p->MapworkSize, sizeof (Int), gb_p->Iwork, Common) ;
    }
  }

  PRINTF("\n\n\nend t_cholmod_super_numeric..\n\n\n");   

  return (Common->status >= CHOLMOD_OK) ;

}



#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
