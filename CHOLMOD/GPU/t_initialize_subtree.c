/* ========================================================================== */
/* === GPU/t_initialize_subtree.c =========================================== */
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
 *   t_initialize_subtree
 *
 * Description:
 *   Contains functions for initializing
 *   subtrees of the elimination tree.
 *
 */


/* includes */
#include "cholmod_template.h"
#include <string.h>
#include <time.h>







/*
 *  Function:
 *    query_gpu
 *
 *  Description:
 *    Queries GPU properties (clock speed, # SMs, etc.) 
 */
void TEMPLATE2(CHOLMOD(query_gpu)) (int *clockRate, int *sm, int *ipc, int gpuid)
{  
#ifdef SUITESPARSE_CUDA
    struct cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, gpuid);
    *clockRate = prop.clockRate; 
    *sm = prop.multiProcessorCount;
    *ipc = 64*2;  /* 64DP ALUs x 2DP (1DP FMA per cycle) */

    PRINTF("GPU Info:\n");
    PRINTFV("\tclock rate:      %d\n",*clockRate);
    PRINTFV("\t# SMs:           %d\n",*sm);
    PRINTFV("\tipc:             %d\n",*ipc);
    PRINTFV("\tname:		%s\n",prop.name);
#endif
}


/*
 *  Function:
 *    query_cpu
 *
 *  Description:
 *    Queries CPU properties (clock speed, # SMs, etc.)
 */
void TEMPLATE2(CHOLMOD(query_cpu)) (int *clockRate, int *sm, int *ipc, int numThreads)
{
    *clockRate = 2300000;	/* clock speed in kilohertz*/
    *sm = numThreads/2;     /* # cores (16 cores) */ 
    *ipc = 16;  		/* 16DP instructions per cycle */

    if(*sm == 0) *sm = 1;

    PRINTF("CPU Info:\n");
    PRINTFV("\tclock rate:	%d\n",*clockRate);
    PRINTFV("\t# cores:		%d\n",*sm);    
    PRINTFV("\tipc:		%d\n",*ipc);
}





/*
 * Function:
 *   binarysearch_tree
 *
 * Description:
 *   Performs binary search to find ideal subtree sizes for
 *   the elimination tree. Splits elimination tree into
 *   subtrees that can be factorized concurrently.
 *
 */
void TEMPLATE2 (CHOLMOD (binarysearch_tree))
  (
    cholmod_common *Common,
    cholmod_sparse *A,
    cholmod_factor *L,
    cholmod_global_pointers *gb_p,
    cholmod_cpu_pointers *cpu_p,
    cholmod_tree_pointers *tree_p,
    Int *LpxSub
  )
{
  /* local variables */
  Int n, nls, numSuper, subtree, subtreeSize, subtreeSizeDiff, subtreeSizePrev, max_factor_size;
  Int *supernode_children_num, *supernode_children_num2, *supernode_children_count2, 
      *supernode_levels_subtree_ptrs, *supernode_subtree_ptrs, *level_num_desc_ptrs,
      *level_descendants_ptrs, *Lpi;
  size_t gpu_memtot, cpu_memtot, size_A;
  int search, binarySearch, counts[3];





  /*
   * Set variables & pointers
   */
  /* set host variables */
  n 		= L->n;
  search	= 0;
  gpu_memtot	= 0;
  cpu_memtot	= 0;

  /* set host pointers */
  Lpi	= cpu_p->Lpi;

  /* set tree pointers */  
  supernode_children_num	= tree_p->supernode_children_num;
  supernode_children_num2	= tree_p->supernode_children_num2;
  supernode_children_count2	= tree_p->supernode_children_count2;
  supernode_levels_subtree_ptrs	= tree_p->supernode_levels_subtree_ptrs;
  supernode_subtree_ptrs		= tree_p->supernode_subtree_ptrs;
  level_num_desc_ptrs		= tree_p->level_num_desc_ptrs;
  level_descendants_ptrs	= tree_p->level_descendants_ptrs;  




  /*
   * Determine whether to use root only:
   * Calculate size of Ai, Ap, Ax. If size is larger
   * than Common->dev_mempool_size (device memory) then
   * use only root tree, not our GPU subtrees.
   */
  nls            = Lpi[L->nsuper] - Lpi[0];
  size_A = (nls + A->ncol + A->nzmax + 2)*sizeof(Int) + A->nzmax*sizeof(double);

  if(size_A >= Common->dev_mempool_size && gb_p->runType != 1) {  
    gb_p->runType = 3;          				/* use only root */    
  }    



	

  /* set maximum BRANCH_SIZE (cutoff) */
  if(gb_p->runType == 1)  subtreeSize = L->xsize + 1;  
  else 			  subtreeSize = L->xsize - 1;           /* initial subtree size (at least one supernode on root alg.) */
  subtreeSizePrev = 0;
  subtreeSizeDiff = 0;
  search = 0;
  subtree = 0;
  binarySearch = (int)(BINARY_SEARCH);

  /* case if factor (subtree size) is larger than GPU memory available */
  if(subtreeSize > (Int)(Common->dev_mempool_size/8)) {
    subtreeSize = (Int)(Common->dev_mempool_size/8);
  }

  



  /* Binary Search loop
   *   conditions:
   *     1. BINARY_SEARCH steps reached
   *     2. factor fits GPU memory
   *     3. factor fits CPU (pinned) memory
   */
  while(search <= binarySearch || gpu_memtot > Common->dev_mempool_size || cpu_memtot > Common->dev_mempool_size) {


    /* case binary search could not find small enough subtree to fit in GPU, use root only.. */
    if ( subtreeSize == 1 ) {
    
      PRINTF("subtreeSize = 1, use root only..\n");

      /* no subtree fits GPU, so use root only */
      gb_p->runType = 3;      

      /* set maximum BRANCH_SIZE (cutoff) */
      subtreeSize = L->xsize - 1;           /* initial subtree size (at least one supernode on root alg.) */
      subtreeSizePrev = 0;
      subtreeSizeDiff = 0;
      search = 0;
      subtree = 0;
    
      /* case if factor (subtree size) is larger than GPU memory available */
      if(subtreeSize > (Int)(Common->dev_mempool_size/8)) {
        subtreeSize = (Int)(Common->dev_mempool_size/8);
      }

    }



    /* clear local variables */
    gb_p->numSubtree  = 0;          /* # subtrees in tree */
    max_factor_size  = 0;          /* max factor size in any subtree */
    gb_p->maxCsize   = 0;          /* max C size in batch & streams in any level in any subtree */
    gb_p->maxndesc   = 0;          /* max # descendants in batch & streams in any level in any subtree */
    gb_p->maxbatch   = 0;	   /* max batch size (# supernodes) in any level */
    counts[0]        = 0;
    counts[1]        = 0;
    counts[2]        = 0;




    /*
     * Store subtrees of tree
     *
     * Description:
     * traverse through tree in two manners:
     * 1. DESCEND if supernode has children that haven't been touched (added to subtree)
     * 2. ASCEND if supernode has no children or if all have been touched (added to subtree)
     *
     * Lastly, store supernodes in last subtree of tree that does not fit GPU
     */
  
    /* copy # children per supernode */
    memcpy(supernode_children_num, supernode_children_num2, L->nsuper*sizeof(Int));
  
    /* reset counters */
    memset(supernode_children_count2, 0, L->nsuper*sizeof(Int)); 
  
    TEMPLATE2 (CHOLMOD (build_subtree))( L,
					gb_p,
					tree_p,
					subtreeSize);

  
  
  
  
    /*
     * pre-processing subtrees
     *
     * Description:
     * 1. Order subtrees by levels
     * 2. Find amount of GPU memory needed for batching
     * 3. Find batching cutoff (up to what level to batch over supernodes)
     */
  
    /* loop over subtrees */
    for(subtree = 0; subtree < gb_p->numSubtree; subtree++) {
  
  
      numSuper = supernode_subtree_ptrs[subtree+1] - supernode_subtree_ptrs[subtree];  /* get # of supernodes in subtree */
  
      /* copy # childrens per supernode */
      memcpy(supernode_children_num, supernode_children_num2, L->nsuper*sizeof(Int));
  
      level_num_desc_ptrs[subtree]          = counts[0];  
      level_descendants_ptrs[subtree]       = counts[1];  
      supernode_levels_subtree_ptrs[subtree] = counts[2];  
  
      /* get # children in root supernode of root tree */
      TEMPLATE2 (CHOLMOD (get_children_root))( Common,
					       gb_p,
					       tree_p,
					       numSuper,
					       subtree);

  
      /* get size of factor (Lx) in current subtree */
      TEMPLATE2 (CHOLMOD (get_factor_size))( gb_p,
					     cpu_p,
					     tree_p,
					     numSuper,
					     subtree,
					     &max_factor_size,
                                             LpxSub);

  
      /* get current subtree size/info */
      TEMPLATE2 (CHOLMOD (process_subtree))( Common,
					    A,
					    L,
					    gb_p,
					    cpu_p,
					    tree_p,
					    n,
					    numSuper,
					    subtree,
					    max_factor_size,
                                            counts);
  
  
    } /* end loop over subtrees */
  
  
  
    /*
     * find amount GPU memory needed for subtreeing
     *
     * Description:
     * calculates total amount of GPU memory needed for current BRANCH_SIZE.
     * If larger then reduce subtree size, if smaller increase subtree size.
     */
    nls             	= Lpi[L->nsuper] - Lpi[0];
    gb_p->LxSize    	= (max_factor_size)*sizeof(double);          	/* size of factor */
    gb_p->CSize     	= (gb_p->maxCsize)*sizeof(double);          	/* size of C buffer */
    gb_p->LsSize    	= (nls+1)*sizeof(Int);                       
    gb_p->MapSize       = (n+1)*sizeof(Int)*(gb_p->maxbatch);          /* size of Map */
    gb_p->ApSize        = (A->ncol+1)*sizeof(Int);                   
    gb_p->AiSize        = A->nzmax*sizeof(Int);                      
    gb_p->AxSize        = A->nzmax*sizeof(double);                   
    gb_p->dimDescSize   = (gb_p->maxndesc)*sizeof(int);               	/* size of dimension arrays for desc. */
    gb_p->ptrDescSize   = (gb_p->maxndesc)*sizeof(double *);          	/* size of pointer arrays for desc. */
    gb_p->dimSuperSize  = sizeof(int)*(gb_p->maxbatch);       		/* size of dimension arrays for super. */
    gb_p->ptrSuperSize  = sizeof(double *)*(gb_p->maxbatch);     	/* size of pointer arrays for super. */
  
  
    /* size of Ap, Ai, Ax buffers (0 if GPU subtrees not used) */
    if(gb_p->runType != 1 && gb_p->runType != 3) size_A = gb_p->ApSize + gb_p->AiSize + gb_p->AxSize;    	/* if not root and not CPU only */
    else	  				 size_A = 0;


    /* total amount of GPU memory needed */
    gpu_memtot = gb_p->LxSize + gb_p->CSize + gb_p->LsSize + gb_p->MapSize + size_A + (14*(gb_p->dimDescSize) + 6*(gb_p->ptrDescSize) + 13*(gb_p->dimSuperSize) + 3*(gb_p->ptrSuperSize)) +
                 2*(gb_p->maxbatch)*sizeof(int) + sizeof(int);
  
    /* total amount of CPU memory needed (pinned memory) */
    cpu_memtot = gb_p->LxSize + (14*(gb_p->dimDescSize) + 6*(gb_p->ptrDescSize) + 13*(gb_p->dimSuperSize) + 3*(gb_p->ptrSuperSize));
  
    /* print memory info */
    PRINTFV("binary step: %d\n",search);
    PRINTFV("\trunType:	 	  %ld \n", gb_p->runType);
    PRINTFV("\tA->nzmax:          %ld \n", A->nzmax);
    PRINTFV("\tLxSize:            %ld \n", gb_p->LxSize);
    PRINTFV("\tCSize:             %ld \n", gb_p->CSize);
    PRINTFV("\tMapSize:           %ld \n", gb_p->MapSize);
    PRINTFV("\tLsSize:            %ld \n", gb_p->LsSize);
    PRINTFV("\tApSize:            %ld \n", gb_p->ApSize);
    PRINTFV("\tAiSize:            %ld \n", gb_p->AiSize);
    PRINTFV("\tAxSize:            %ld \n", gb_p->AxSize);
    PRINTFV("\tbatch lists:       %ld \n", (14*(gb_p->dimDescSize) + 6*(gb_p->ptrDescSize) + 13*(gb_p->dimSuperSize) + 3*(gb_p->ptrSuperSize)) + 2*(gb_p->maxbatch)*sizeof(int));
    PRINTFV("\tcpu_mem_available: %ld \n",Common->dev_mempool_size);
    PRINTFV("\tcpu_mem_used:      %ld \n",cpu_memtot);
    PRINTFV("\tgpu_mem_available: %ld \n",Common->dev_mempool_size);
    PRINTFV("\tgpu_mem_used:      %ld \n",gpu_memtot);
    PRINTF("\n\n");
  
  
  
    /*
     * Update size of subtree.
     *
     * Update subtreeSize by subtreeSizeDiff amount, where subtreeSizeDiff is half the difference
     * between the previous and current size. Increase if the subtreeSize is smaller than the available
     * GPU (or CPU) memory, and decrease otherwise. Also store the current subtree size as subtreeSizePrev.
     */
  
    /* Subtree size change to update. Half the difference between the previous and current subtree size. */
    subtreeSizeDiff = (Int)((float)(labs(subtreeSize - subtreeSizePrev))/2.0 + 0.5);
  
    /* Do not let it exceed half the subtree size. */
    if ( subtreeSizeDiff > (subtreeSize)/2) {
      subtreeSizeDiff = (subtreeSize)/2 ;
    }
  
    /* store previous subtree size */
    subtreeSizePrev = subtreeSize;
  
    /* update size of subtree */
    /* case if exceed GPU or CPU memory, reduce subtree size */
    if (gpu_memtot > Common->dev_mempool_size || cpu_memtot > Common->dev_mempool_size) {
      subtreeSize -= subtreeSizeDiff;
    }
    /* case if not exceed GPU nor CPU memory, increase subtree size */
    else {
      subtreeSize += subtreeSizeDiff;
    }
  
  
    /* break conditions for exiting binary search loop:
     *    1. BINARY_SEARCH steps reached
     *    2. GPU mem does not exceed limit
     *    3. if subtree size converges (BRANCH_SIZE_DIFF == 0)
     *    4. if subtree size reaches size of factor or size of factor - 1 (depends on SINGLE_BRANCH)
     *    5. if root_only, subtree defaults to only root algorithm
     *    6. if CPU_only
     */
    if(((gpu_memtot < Common->dev_mempool_size && cpu_memtot < Common->dev_mempool_size) &&
       (search >= binarySearch || !(subtreeSizeDiff) || subtreeSize >= L->xsize)) || (gb_p->runType == 1) || (gb_p->runType == 3))       
         {
           break;
         }
  

    /* increment binary step */
    search++;

  } /* end binary search loop*/


}










/*
 * Function:
 *   loadbalance_gpu
 *
 * Description:
 *   Load balances subtrees on multiple devices. Four cases:
 *   1. CPU only: sends all subtrees to CPU device (with id CHOLMOD_DEVICE_GPU)
 *   2. root only: sends all subtrees to root (with id CHOLMOD_DEVICE_GPU+1)
 *   3. GPU only: sends subtrees to GPU device (id from 0 to CHOLMOD_DEVICE_GPU-1) & root (id CHOLMOD_DEVICE_GPU+1) 
 *   4. hybrid:   sends subtrees to GPU device (id from 0 to CHOLMOD_DEVICE_GPU-1), CPU device (id CHOLMOD_DEVICE_GPU) & root (id CHOLMOD_DEVICE_GPU+1)
 *
 *   The load-balancing algorithm computes the total work on each device (runtime of all subtrees computed as flop/flops). Then it assigns each subtree
 *   to the device with least amount of work in a cyclic fashion.
 */  
void TEMPLATE2 (CHOLMOD (loadbalance_gpu))
  (
    cholmod_common *Common,
    cholmod_global_pointers *gb_p,
    cholmod_tree_pointers *tree_p,
    cholmod_loadbalance_pointers *lb_p
  )
{
  /* structure for device properties */
  typedef struct props{
    int clockRate;
    int sm;
    int ipc;
  };

  /* local variables */
  double GPUflops, CPUflops, flop, GPUtime, CPUtime;
  double *subtreeSize, *supernode_flop, *workPerDevice;
  int i, j, runType, numDevice, numSubtree;
  Int s;
  Int *supernode_subtree, *supernode_subtree_ptrs, *numSubtreePerDevice, *listSubtreePerDevice;
  struct props gpu, cpu; 
  struct cholmod_subtree_order_t *subtreeReorder;





  /*
   * Set variables & pointers
   */
  /* set variables */
  runType	= gb_p->runType;  
  numSubtree 	= gb_p->numSubtree;  

  /* set load-balance pointers */
  subtreeSize		= lb_p->subtreeSize;
  numSubtreePerDevice	= lb_p->numSubtreePerDevice;
  listSubtreePerDevice	= lb_p->listSubtreePerDevice;
  subtreeReorder		= lb_p->subtreeReorder;
  workPerDevice		= lb_p->workPerDevice;

  /* set tree */
  supernode_flop	= tree_p->supernode_flop;
  supernode_subtree	= tree_p->supernode_subtree;
  supernode_subtree_ptrs	= tree_p->supernode_subtree_ptrs;





  /* get number of devices to use:
   *   1. GPU only: Common->numGPU
   *   2. hybrid:   Common->numGPU+1  
   */
  if(runType == 1)		numDevice = 1;				/* CPU only */
  else if(runType == 2)		numDevice = Common->numGPU;		/* GPU only */
  else				numDevice = Common->numGPU + 1;		/* GPU + CPU (hybrid) */





  /* compute theoretical GPU & CPU flops */
  TEMPLATE2(CHOLMOD(query_gpu)) (&(gpu.clockRate), &(gpu.sm), &(gpu.ipc), 0);
  TEMPLATE2(CHOLMOD(query_cpu)) (&(cpu.clockRate), &(cpu.sm), &(cpu.ipc), Common->ompNumThreads);
  GPUflops = (double)(gpu.clockRate*gpu.sm*gpu.ipc)/(double)(1.0e+6);		/* GPU peak theoretical performance (in gflops) */
  CPUflops = (double)(cpu.clockRate*cpu.sm*cpu.ipc)/(double)(1.0e+6);           /* CPU peak theoretical performance (in gflops) */
  PRINTFV("GPU peak flops rate: %f\n",GPUflops);
  PRINTFV("CPU peak flops rate: %f\n",CPUflops); 





  /* Store subtree info (size and id): 
   * computes the cumulative number of floating-point operations (flop) in each subtree. 
   */
  for(i = 0; i < numSubtree; i++) 
  {

    subtreeReorder[i].id = i;                            /* subtree id */
    int numSuper = supernode_subtree_ptrs[i+1] - supernode_subtree_ptrs[i]; 

    /* loop over supernodes in subtree */
    for(j = 0; j < numSuper; j++) {
      s = supernode_subtree[supernode_subtree_ptrs[i] + j];
      subtreeReorder[i].size += supernode_flop[s];	/* subtree size (# flop in all its supernodes) */   
    }

    /* convert to gflop */
    subtreeReorder[i].size *= 1.0e-9;	
    subtreeSize[i] = subtreeReorder[i].size;
  }





  /* reorder subtrees by size (largest to smallest number of flop) */
  qsort(subtreeReorder, numSubtree, sizeof(struct cholmod_subtree_order_t),CHOLMOD(subtree_comp));





  /* Set subtrees for each device
   * Finds the device with least work and then adds current subtree to it.
   * The amount of work in the device (workPerDevice) is computed as the total runtime (flop/flop rate)
   * of all the subtrees in the device. The flop rate is the theoretical peak flops of the device (GPU or CPU).  
   */
  PRINTF("\nSubtree Info:\n");
  /* loop over subtrees */
  for(i = 0; i < numSubtree; i++) 
  {

    int minDevice = 0;
    double min, size;

    /* set initial device */
    min = workPerDevice[0];


    /* case CPU device (CPU only) */
    if(runType == 1)    
    {
      minDevice = Common->numGPU;                          	/* set CPU device */
    }
    /* case root (last subtree) */
    else if(subtreeReorder[i].id == numSubtree-1) 
    {
      minDevice = Common->numGPU + 1;				/* set root */
    } 
    /* case GPU or CPU device (GPU only or hybrid) */
    else 
    {
      /* find device with least work */      
      for(j = 1; j < numDevice; j++) {
        if(min > workPerDevice[j]) {
          min = workPerDevice[j];
          minDevice = j;					/* set GPU or CPU device */
        }
      }
    }


    /* compute size (execution time) for subtree */
    flop = subtreeReorder[i].size;				/* floating-point operations in subtree */
    GPUtime = flop/GPUflops;					/* GPU runtime */
    CPUtime = flop/CPUflops;					/* CPU runtime */
    if(minDevice == Common->numGPU) 	size = CPUtime;
    else			      	size = GPUtime; 


    /* print subtree info */
    PRINTFV("device:%d ",minDevice);
    PRINTFV("subtree:%d ",subtreeReorder[i].id);
    PRINTFV("workPerDevice:%f ",workPerDevice[minDevice]);
    PRINTFV("subtreeSize:%f ",subtreeReorder[i].size);
    PRINTFV("GPU time:%f ",GPUtime);
    PRINTFV("CPU time:%f\n",CPUtime);


    /* set subtree for selected device (GPU,CPU,root) */
    listSubtreePerDevice[(numSubtreePerDevice[minDevice]++) + minDevice*numSubtree] = (Int)(subtreeReorder[i].id);
    workPerDevice[minDevice] += size;

  } /* end loop over subtrees */





  /* issue less GPUs if not sufficient subtrees */
  if(numSubtree-1 < Common->numGPU) 
  {
    gb_p->numGPU = numSubtree-1;
  }
  else 
  {
    gb_p->numGPU = Common->numGPU;
  }

}










/*
 * Function:
 *   initialize_gpu
 *
 * Description:
 *   initializes for GPU algorithm. 
 *
 */
void TEMPLATE2 (CHOLMOD (initialize_gpu))
  (
    cholmod_common *Common,
    cholmod_factor *L,
    cholmod_sparse *A,
    cholmod_global_pointers *gb_p,
    cholmod_gpu_pointers *gpu_p,
    cholmod_cpu_pointers *cpu_p
  )
{
#ifdef SUITESPARSE_CUDA
   int i, runType, numGPU;
   Int s;


  /* set variables */
  runType	= gb_p->runType;  
  numGPU	= gb_p->numGPU;




  /* initialize GPU (set pointers, copy memory, etc.) 
   * only if there are GPU subtrees  */
  if(runType != 1 && runType != 3) {  
    #pragma omp parallel num_threads(numGPU)
    {
      /* get GPU id (omp thread id) */
      int gpuid = omp_get_thread_num();

      /* set GPU device */
      cudaSetDevice(gpuid);

      /* initialize GPU (set pointers, copy memory, etc.) */
      TEMPLATE2 (CHOLMOD (gpu_init))( Common,
  				      L,
				      A,
				      gb_p,
				      gpu_p,
				      gpuid);
    }
  }

  /* Ensure that all GPU initializations are complete */
  cudaDeviceSynchronize();

#endif
}










/*
 * Function:
 *   initialize_cpu
 *
 * Description:
 *   initializes for root and CPU algorithm.
 *
 */
void TEMPLATE2 (CHOLMOD (initialize_cpu))
  (
    cholmod_common *Common,
    cholmod_factor *L,
    cholmod_global_pointers *gb_p,
    cholmod_cpu_pointers *cpu_p,
    cholmod_tree_pointers *tree_p
  )
{
   int i, runType, numSubtree, numThreads;
   Int s;
   Int *supernode_subtree, *supernode_subtree_ptrs;
   size_t CSize, MapSize;


  /* set variables */
  runType	= gb_p->runType;  
  numSubtree     = gb_p->numSubtree;
  numThreads	= Common->ompNumThreads;
  CSize         = (gb_p->CSize);
  MapSize       = (gb_p->MapSize);

  /* set tree pointers */
  supernode_subtree      = tree_p->supernode_subtree;
  supernode_subtree_ptrs = tree_p->supernode_subtree_ptrs;




  /* set size of Cbuff */
  if(CSize < Common->numGPU*Common->devBuffSize)
    CSize = (Common->numGPU+1)*Common->devBuffSize;

  if(MapSize < (size_t)(Common->numGPU*L->n*sizeof(Int)))
    MapSize = (size_t)(Common->numGPU*L->n*sizeof(Int));




  /* clear Lx factor (supernodes used for root alg.) */
  Int *lpx = L->px;
  #pragma omp parallel for num_threads(numThreads) private(i)
  for(i=supernode_subtree_ptrs[numSubtree-1]; i<supernode_subtree_ptrs[numSubtree]; i++) {
    s = supernode_subtree[i];
    double *ps = (double *)&cpu_p->Lx[lpx[s]];
    memset(ps, 0, sizeof(double));
  }


  /* allocate memory for Cbuff & Map (for factorize_cpu_parallel) */
/*  if (runType != 3 && runType != 2) */
  {

    gb_p->CworkSize = (CSize + sizeof(double) - 1)/sizeof(double);
    gb_p->MapworkSize = (MapSize + sizeof(Int) - 1)/sizeof(Int);

    /* allocate workspace */
    gb_p->Cwork     = CHOLMOD(malloc) (gb_p->CworkSize, sizeof (double), Common) ;
    gb_p->Mapwork   = CHOLMOD(malloc) (gb_p->MapworkSize, sizeof (Int), Common) ;


    /* check if enough memory */
    if (Common->status < CHOLMOD_OK)
    {
      gb_p->Iwork = CHOLMOD(free) (gb_p->IworkSize, sizeof (Int), gb_p->Iwork, Common) ;
      gb_p->Xwork = CHOLMOD(free) (gb_p->XworkSize, sizeof (double), gb_p->Xwork, Common) ;
      gb_p->Bwork = CHOLMOD(free) (gb_p->BworkSize, sizeof (struct cholmod_subtree_order_t), gb_p->Bwork, Common) ;
      gb_p->Cwork = CHOLMOD(free) (gb_p->CworkSize, sizeof (double), gb_p->Iwork, Common) ;
      gb_p->Mapwork = CHOLMOD(free) (gb_p->MapworkSize, sizeof (Int), gb_p->Iwork, Common) ;
      return (FALSE) ;
    }

    cpu_p->C      = gb_p->Cwork;
    cpu_p->Map    = gb_p->Mapwork;

  }


}










/*
 *  Function:
 *    gpu_copy_supernode
 *  
 *  Description:
 *    builds initial elimination tree 
 *
 */
void TEMPLATE2 (CHOLMOD (build_tree))
  (
   cholmod_common *Common,
   cholmod_factor *L,
   cholmod_global_pointers *gb_p, 
   cholmod_cpu_pointers *cpu_p,
   cholmod_tree_pointers *tree_p
   )
{
  /* local variables */
  Int s, k1, k2, nscol, nsrow, psi, psend, ndcol, ndrow, ndrow1, ndrow2, pdx1, pdi1, 
      d, dlarge, kd1, kd2, pdi, pdend, pdi2, dancestor, sparent, id, totdesc, idescendant;
  Int *Super, *SuperMap, *Lpi, *Ls, *Head, *Next, *Lpos, *supernode_root, *supernode_children,
      *supernode_children_count, *supernode_children_num, *supernode_children_ptrs,
      *supernode_parent, *supernode_size, *supernode_size_desc, *ndescendants;
  int childrenPtrs = 0, count;
  double syrkflops,gemmflops,potrfflops,trsmflops;
  double *supernode_flop;




  /* set host pointers */
  Super		= cpu_p->Super;
  SuperMap	= cpu_p->SuperMap;
  Lpi		= cpu_p->Lpi;
  Ls		= cpu_p->Ls;
  Head		= cpu_p->Head;
  Next		= cpu_p->Next; 
  Lpos		= cpu_p->Lpos;

  /* set tree pointers */
  supernode_root		= tree_p->supernode_root;
  supernode_children		= tree_p->supernode_children;
  supernode_children_num	= tree_p->supernode_children_num;
  supernode_children_count	= tree_p->supernode_children_count;
  supernode_children_ptrs	= tree_p->supernode_children_ptrs;
  supernode_parent		= tree_p->supernode_parent;
  supernode_size		= tree_p->supernode_size;
  supernode_size_desc		= tree_p->supernode_size_desc;
  ndescendants			= tree_p->ndescendants;
  supernode_flop		= tree_p->supernode_flop;




  /*
   * Get info of tree:
   * Gathers info of the tree. Visit all
   * supernoeds and collects three things:
   *   1. size
   *   2. parent
   *   3. # children  
   *
   */
  /* loop over supernodes */
  for(s = 0; s < L->nsuper; s++) {  

    /* clear variables */
    id=0;
    totdesc=0;
    idescendant=0;
    syrkflops = 0.0;
    gemmflops = 0.0;
    potrfflops = 0.0;
    trsmflops = 0.0;

    /* get supernode dimensions */
    k1 = Super[s] ;
    k2 = Super[s+1] ;
    nscol = k2 - k1 ;
    psi = Lpi[s] ;
    psend = Lpi[s+1];
    nsrow = psend - psi;

    /* store maximum nsrow & nscol in tree*/
    if(nsrow > gb_p->maxnsrow) gb_p->maxnsrow = nsrow;
    if(nscol > gb_p->maxnscol) gb_p->maxnscol = nscol;

    /* get number of descendants in supernode */
    TEMPLATE2 (CHOLMOD (gpu_num_descendants))( Common,
					       cpu_p,
					       tree_p,
                                               s);

    /* get current supernode */
    d = Head[s];
    dlarge = Next[d];



    /* loop over descendants of supernode */
    while(idescendant < ndescendants[s]) {

      /* get current descendant */
      if (idescendant > 0) {
        d = dlarge;
        dlarge = Next[dlarge];
      }

      /* increment descendant count */
      idescendant++;

      /* get descendant dimensions */
      kd1 = Super [d] ;
      kd2 = Super [d+1] ;
      ndcol = kd2 - kd1 ;
      pdi = Lpi [d] ;
      pdend = Lpi [d+1] ;
      ndrow = pdend - pdi ;
      pdi1 = pdi + Lpos[d];
      for (pdi2 = pdi1 ; pdi2 < pdend && Ls [pdi2] < k2 ; pdi2++) ;
      ndrow1 = pdi2 - pdi1 ;
      ndrow2 = pdend - pdi1 ;

      /* get next descendant */
      Lpos [d] = pdi2 - pdi ;
      if (Lpos [d] < ndrow) {
        dancestor = SuperMap [Ls [pdi2]] ;
        Next [d] = Head [dancestor] ;
        Head [dancestor] = d ;
      }

      /* cumulate total size of all descendants in current supernode */
      totdesc += ndrow2*ndrow1;

      /* compute syrk & gemm flops in current descendant */
      syrkflops += (double)(ndrow1*ndrow1*ndcol);
      gemmflops += (double)(2.0*(ndrow2-ndrow1)*ndrow1*ndcol);

      id++;

    } /* end loop over descendants */

    /* compute potrf & trsm flops in current supernode */
    potrfflops = (double)(nscol*nscol*nscol/3.0);
    trsmflops = (double)((nsrow-nscol)*nscol*nscol);


    /* get next supernode */
    if(nsrow > nscol) {
      Lpos [s] = nscol ;
      sparent = SuperMap [Ls [psi + nscol]] ;
      Next [s] = Head [sparent] ;
      Head [sparent] = s ;
    }
    Head [s] = EMPTY ;


    /* store tree information */
    supernode_size_desc[s] = totdesc;   	    		/* store total size of all descendants in current supernode */
    supernode_size[s] += totdesc;   				/* store total size of current supernode */
    supernode_flop[s] = syrkflops+gemmflops+potrfflops+trsmflops;	/* store total flops in current supernode */
    if(nsrow > nscol) {		/* case if supernode has parent */
      sparent = SuperMap[Ls [psi + nscol]] ;        
      supernode_size[sparent] += supernode_size[s]; 	/* add supernode's size to its parent */
      supernode_parent[s] = sparent;                		/* store supernode's parent */
      supernode_children_num[sparent]++;            		/* increment # of children of supernode's parent */
    }
    else {			/* case if supernode has no parent */
      supernode_parent[s] = EMPTY;
    }

  } /* end loop over supernodes */





  /*
   * Store children of tree:
   * Builds elimination tree. Visits
   * all supernoes and stores their
   * children. 
   */
  /* loop over supernodes */
  for(s = 0; s < L->nsuper; s++) {

    sparent = supernode_parent[s]; 

    if(sparent > 0) {		/* case if supernode has parent */

      count = supernode_children_count[sparent];         /* get # children the supernode's parent has */

      if(!count) {		/* case if supernode does not have siblings (its parent has no other descendants) */
        supernode_children_ptrs[sparent] = childrenPtrs; /* set children pointer to child */
      }

      /* store children info */
      id = supernode_children_ptrs[sparent] + count;     /* index to store child (# siblings) */
      supernode_children[id] = s;                        /* store supernode as a child */
      supernode_children_count[sparent]++;               /* increment # siblings of supernode (or # children (descendants) parent has) */

      if(!count)
        childrenPtrs += supernode_children_num[sparent]; /* increment pointer to children */
    }
    else {			/* case if supernode has no parent (it is the root of a tree) */
      supernode_root[(gb_p->numRoot)++] = s;                  /* store roots of trees */
    }

  } /* end loop over supernodes */


}










/*
 *  Function:
 *    build_subtree
 *
 *  Description:
 *    builds a subtree of the elimination tree
 *
 */
void TEMPLATE2 (CHOLMOD (build_subtree))
  (
   cholmod_factor *L,
   cholmod_global_pointers *gb_p,
   cholmod_tree_pointers *tree_p,
   Int subtreeSize
   )
{
  /* local variables */
  Int i, j, s, id;
  Int *supernode_root, *supernode_children, *supernode_children_ptrs, *supernode_children_num, *supernode_children_count2,
      *supernode_parent, *supernode_subtree, *supernode_subtree_ptrs, *supernode_size;
  int subtree, first, numRoot, runType;
 
  /* set variables */
  j 			= 0; 
  gb_p->numSubtree 	= 0;
  numRoot 		= gb_p->numRoot;
  runType		= gb_p->runType;  

  /* set tree pointers */
  supernode_root		= tree_p->supernode_root;
  supernode_children		= tree_p->supernode_children;
  supernode_children_ptrs	= tree_p->supernode_children_ptrs;
  supernode_children_num	= tree_p->supernode_children_num;
  supernode_children_count2	= tree_p->supernode_children_count2;
  supernode_parent		= tree_p->supernode_parent;
  supernode_subtree		= tree_p->supernode_subtree;
  supernode_subtree_ptrs		= tree_p->supernode_subtree_ptrs; 
  supernode_size                = tree_p->supernode_size;


  /* 
   * Build subtrees of tree:
   *
   * traverse tree and store supernodes in 
   * corresponding subtrees. Use depth-first
   * traversal.
   * Steps:  
   *   1. start from root supernode (there can be
   *      multiple roots)
   *
   *   2. descend tree until base (leaves) of tree reached
   *      (where supernodes have no children)
   *
   *   3. add supernode to subtree (subtree) if:
   *        a. supernode size < subtree size
   *        b. supernode has no children, or all
   *           children have already been visited
   *
   *   4. ascend tree if: 
   *        a. supernode has no children
   *        b. all children visited
   *
   * Use variable 'first' to determine when a
   * new subtree starts. Set 'first' to head 
   * (root) of subtree (subtree) and stop when 
   * 's' =='first', which means we've returned
   * to the head supernode of the subtree. 
   * Increment the # subtrees whenever this happens. 
   * */

  /* loop over all roots (trees) */
  for(i=0; i<numRoot; i++) {

    s = tree_p->supernode_root[i];      /* set root of tree */
    first = 0;



    /* loop: traverse (depth-first) through supernodes of current tree */
    while(1) {

      /* if returned to first supernode in subtree (exit condition) */
      if(first == s) {
        first = 0;
      }


      /* case: store supernodes in subtree only if: 
       *       1. size of supernode is smaller than size of subtree
       *          note that supernode_size contains size of all its descendants (children)
       *       2. there are at least SUPERNODE_MIN supernods in the elimination tree (root_only = 0)
       *       3. Ai,Ap,Ax are small enough to fit device (root_only = 0)
       */
      if(supernode_size[s] <= subtreeSize && (runType != 3) && (runType != 1)) {      

        /* case: if first supernode in subtree (root of subtree) */
        if(!first) {
          first = supernode_parent[s];                	/* store first supernode */
          supernode_subtree_ptrs[(gb_p->numSubtree)++] = j;  	/* set pointer to current subtree */
        }

        /* case: if supernode has no children */
        if(supernode_children_count2[s]==0) {
          supernode_subtree[j++] = s;                  	/* store supernode into subtree */
        }
      }


      /* case: descend to next child (traverse down tree)
       *       if supernode has children that haven't been added to the subtree
       */
      if(supernode_children_count2[s] < supernode_children_num[s]) {
        id = supernode_children_ptrs[s] + supernode_children_count2[s];       /* get id of children of the supernode */
        supernode_children_count2[s]++;                                       /* increment children count */
        s = supernode_children[id];                                           /* get id of the child (descendant) */
      }

      /* case: ascend to parent (traverse up tree)
       *       if supernode has no children or all children have been added to subtree
       *       and supernode is not a root (since roots have no parents..)
       */
      else if (s != supernode_root[i]){
        s = supernode_parent[s];                                         /* get id of the parent */
      }

      /* exit if root of tree reached */
      if(s == supernode_root[i] && supernode_children_count2[s] == supernode_children_num[s]) {
        break;
      }

    } /* end loop to traverse tree */
  } /* end loop over trees */





  /* 
   * Build last (root) subtree of tree:
   *
   * store supernodes that do not fit GPU ( > subtree size)
   * into last subtree (root subtree). These supernodes are
   * typically located at the top of the tree. 
   */

  supernode_subtree_ptrs[(gb_p->numSubtree)++] = j;  	/* set poiner to last subtree */

  for(s=0; s < L->nsuper; s++) {                 		/* loop over supernodes */
    /* case if size of candidate subtree > cutoff subtree size */
    if(supernode_size[s] > subtreeSize || (runType == 3) || (runType == 1)) {    
      supernode_subtree[j++] = s;              	/* store supernode in subtree */
    }
  } /* end loop over supernodes */

  /* set pointer for end of last subtree */
  supernode_subtree_ptrs[gb_p->numSubtree] = j;
 
}










/*
 *  Function:
 *    get_children_root
 *
 *  Description:
 *    stores # children on root supernodes for last (root) subtree
 *    builds a subtree of the elimination tree
 *
 */
void TEMPLATE2 (CHOLMOD (get_children_root))
  (
   cholmod_common *Common,
   cholmod_global_pointers *gb_p,
   cholmod_tree_pointers *tree_p,
   Int numSuper,
   Int subtree
   )
{
  /* local variables */
  Int i, j, k, s;
  Int *supernode_children, *supernode_children_ptrs, *supernode_subtree, *supernode_subtree_ptrs, *supernode_children_num;
  int num, child, numSubtree, numThreads; 

  /* set variables */
  numSubtree 	= gb_p->numSubtree;
  numThreads	= Common->ompNumThreads;

  /* set tree poitners */
  supernode_children		= tree_p->supernode_children;
  supernode_children_ptrs	= tree_p->supernode_children_ptrs;
  supernode_subtree		= tree_p->supernode_subtree;
  supernode_subtree_ptrs		= tree_p->supernode_subtree_ptrs;
  supernode_children_num	= tree_p->supernode_children_num;
  


  /* 
   * Get children on root supernodes:  
   * get # of children on root supernodes for root (last) subtree. 
   *
   */

  /* case if root (last) subtree */
  if(subtree == numSubtree-1 ) {

    /* loop over supernodes */
#pragma omp parallel for num_threads(numThreads)
    for(i = 0; i < numSuper; i++) {

      /* get supernode id */
      s = supernode_subtree[supernode_subtree_ptrs[subtree] + i];
      num = 0;

      /* loop over children of supernode */
      for(j = 0; j < supernode_children_num[s]; j++) {

        /* get children id */
        child = supernode_children[supernode_children_ptrs[s] + j];

        /* loop over supernodes in last subtree */
        for(k=0; k < numSuper; k++) {
          if(child == supernode_subtree[supernode_subtree_ptrs[subtree] + k]) { /* case: is it a child? */
            num++;
          }
        } /* end loop over supernodes */

      } /* end loop over children*/

      supernode_children_num[s] = num;  /* store # children supernode has in last subtree */

    } /* end supernode loop */
  }

}










/*
 *  Function:
 *    process_subtree
 *    
 *  Description:
 *    processes a subtree of the elimination tree. Stores supernodes
 *    in levels, and finds the maximum batch size for each level, 
 *    given a fixed amount of device memory.
 * 
 */
void TEMPLATE2 (CHOLMOD (process_subtree))
  (
   cholmod_common *Common,
   cholmod_sparse *A,
   cholmod_factor *L,
   cholmod_global_pointers *gb_p,
   cholmod_cpu_pointers *cpu_p,
   cholmod_tree_pointers *tree_p,
   Int n,
   Int numSuper,
   Int subtree,
   Int max_factor_size,  
   int *counts)
{
  /* local variables */
  int batchdescflag, desc, count0, count1, count2, nsupernodes, stream, batch, Csize, ndesc, 
      maxsubtreeCsize, maxsubtreendesc, maxsubtreebatch, maxnumdescendantsperlevel, nbatch,
      maxsubtreeCsize_prev, maxsubtreendesc_prev, maxsubtreebatch_prev, maxnumdescendantsperlevel_prev, nbatch_prev, runType;
  Int s, i, processed_nodes, node,  num_levels, sparent;
  Int *Lpi, *supernode_levels, *supernode_levels_ptrs, *supernode_levels_subtree_ptrs, *supernode_subtree, *supernode_subtree_ptrs,
      *level_descendants, *supernode_children_num, *supernode_parent, *supernode_size_desc, *supernode_num_levels,
      *level_num_desc, *supernode_batch, *ndescendants;
  size_t nls, LxSize, CSize, LsSize, MapSize, ApSize, AiSize, AxSize, dimDescSize, ptrDescSize, dimSuperSize, ptrSuperSize,
         gpu_memtot, gpu_memtot_prev;





  /* set variables */
  processed_nodes = 0;
  num_levels 	= 0;
  count0 	= counts[0];
  count1 	= counts[1];
  count2 	= counts[2];
  runType	= gb_p->runType;  

  /* set host pointers */
  Lpi 	= cpu_p->Lpi; 

  /* set tree pointers */
  supernode_levels		= tree_p->supernode_levels;
  supernode_levels_ptrs		= tree_p->supernode_levels_ptrs;
  supernode_levels_subtree_ptrs	= tree_p->supernode_levels_subtree_ptrs;
  supernode_subtree		= tree_p->supernode_subtree;
  supernode_subtree_ptrs		= tree_p->supernode_subtree_ptrs;
  level_descendants		= tree_p->level_descendants;
  supernode_children_num	= tree_p->supernode_children_num;
  supernode_parent		= tree_p->supernode_parent;
  supernode_size_desc		= tree_p->supernode_size_desc;
  supernode_num_levels		= tree_p->supernode_num_levels;
  level_num_desc		= tree_p->level_num_desc;
  supernode_batch		= tree_p->supernode_batch;
  ndescendants			= tree_p->ndescendants;




  /*
   * Process subtree:  
   * First store all supernodes within
   * levels. Then visit all supernodes
   * in a level and get the amount of
   * memory needed for batching them. 
   */ 
 
  /* loop over levels in subtree (until all supernodes are processed) */
  while(processed_nodes != numSuper) {
    
    /* reset variables */
    nsupernodes = 0;
    batchdescflag = 0;
    desc = 0;
    stream = 0;
    batch = 0;
    Csize = 0;
    ndesc = 0;
    maxsubtreeCsize = 0;
    maxsubtreendesc = 0;



 

    /* Store supernods in current level: 
     * This just involves selecting supernodes that have no
     * children (belong to the current level).
     */

    /* pointer to levels in subtree */
    supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+num_levels] = count2;
    nsupernodes = 0;

    /* loop over supernodes */
    for(i=0; i < numSuper; i++) {  

      s = supernode_subtree[supernode_subtree_ptrs[subtree] + i]; 

      /* store supernodes that belong to current level */
      if(supernode_children_num[s] == 0) { /* case supernode has no children (belongs to current level) */

        supernode_levels[count2++] = s;    	/* store supernode in level */                   
        processed_nodes++;                 		/* increment processed supernode coutner */
        nsupernodes++;                     		/* increment # supernodes in level */

      }
    } /* end loop over supernodes */





    /* 
     * Update supernodes: 
     * Remove supernodes in current level
     * from their parent's children list.
     *
    */

   /* loop over supernodes in level */
    for(i = 0; i < nsupernodes; i++) {

      node = supernode_levels[supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+num_levels]+i];	/* get supernode */
      supernode_children_num[node] = EMPTY;   									/* empty children of supernode */
      sparent = supernode_parent[node];   									/* get parent of supernode*/

      /* case if parent of supernode has children */
      if(sparent != EMPTY) {
        supernode_children_num[sparent]--;  									/* remove supernode as child of its parent (supernode is processed) */
      }

      /* get maximum # of descendants in a supernode in current level */
        if(ndescendants[node] > desc) {
          desc = ndescendants[node];
        }

    } /* end loop over supernodes */


    /* store maximum # descendants a supernode has in current level */
    level_descendants[count1++] = desc;





    /*
     * Get batching info:
     *
     * Compute the amount of memory needed for batching
     * supernodes. That is, store three variables:
     *
     *   1. maxbatch: 
     *     a. maximum batch size (of supernodes) in any given level
     *     b. size of buffers to store lists of supernode dimensions  
     *
     *   2. maxndesc:
     *     a. maximum number of descendants in any batch
     *     b. size of buffers to store lists of descendant dimensions
     *
     *   3. maxCsize: 
     *     a. maximum cumulative size of descendants in a batch
     *     b. size of buffer to store schur complements
     *     
     * The algorithm below finds the optimal (largest) batch size (# supernodes)
     * to be used for each level. 
     *
     * But only do this if the subtree is not the root subtree (the last top-of-tree 
     * subtree).
     *
     */


    /* 
     * case if: 
     *   1. one of GPU subtrees (not root subtree)
     *   2. not root only
     *   3. not CPU only 
     */
    if((subtree != gb_p->numSubtree-1) && (runType != 3) && (runType != 1)) {    



      /* reset variables */
      maxsubtreeCsize = (int)gb_p->maxCsize;
      maxsubtreendesc = (int)gb_p->maxndesc;
      maxsubtreebatch = (int)gb_p->maxbatch;
      maxnumdescendantsperlevel = 0;
      gpu_memtot = 0;
      gpu_memtot_prev = gpu_memtot;
      nbatch = 1;
     
      /* while loop to find batch size for current level */
      while(1) {

        /* reset variables */
        Csize = 0;
        ndesc = 0;
        gpu_memtot_prev = gpu_memtot;
        maxsubtreeCsize_prev = maxsubtreeCsize;
        maxsubtreendesc_prev = maxsubtreendesc;
        maxsubtreebatch_prev = maxsubtreebatch;
        maxnumdescendantsperlevel_prev = maxnumdescendantsperlevel;

 
        /* loop over supernodes in level */
        for(i = 0; i < nsupernodes; i++) {

          /* get supernode */
          node = supernode_levels[supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+num_levels]+i];    


          /* reset variables (new batch) */
          if( !(i % nbatch ) ) {
            Csize = 0;
            ndesc = 0;       
          }

          Csize += supernode_size_desc[node];			/* add to total size of C buffer needed for storing schur complement in current batch of supernodes */
          ndesc += ndescendants[node];				/* add to total # descendants in current batch of supernodes */

          if(Csize > maxsubtreeCsize) maxsubtreeCsize = Csize;		/* store maximum C buffer size in any given level */
          if(ndesc > maxsubtreendesc) maxsubtreendesc = ndesc;		/* store maximum # descendants in any given level */
	  if(nbatch > maxsubtreebatch) maxsubtreebatch = nbatch;		/* store maximum batch size in any given level */
          if(ndesc > maxnumdescendantsperlevel) maxnumdescendantsperlevel = ndesc;	


        } /* end loop over supernodes in level */
 
        /* find amount GPU memory needed for subtreeing */
        nls            = Lpi[L->nsuper] - Lpi[0];
        LxSize         = max_factor_size*sizeof(double);        /* size of factor */
        CSize          = maxsubtreeCsize*sizeof(double);         /* size of C buffer */
        LsSize         = (nls+1)*sizeof(Int);                   
        MapSize        = (n+1)*sizeof(Int)*(maxsubtreebatch);    /* size of Map */
        ApSize         = (A->ncol+1)*sizeof(Int);               
        AiSize         = A->nzmax*sizeof(Int);                  
        AxSize         = A->nzmax*sizeof(double);               
        dimDescSize    = maxsubtreendesc*sizeof(int);            /* size of list of dimensions for descendants */
        ptrDescSize    = maxsubtreendesc*sizeof(double *);       /* size of list of pointers for descendants */
        dimSuperSize   = sizeof(int)*(maxsubtreebatch);          /* size of list of dimensions for supernodes */
        ptrSuperSize   = sizeof(double *)*(maxsubtreebatch);     /* size of list of pointers for supernodes */

        /* compute total amount of GPU memory needed */
        gpu_memtot_prev = gpu_memtot;
        gpu_memtot = LxSize + CSize + LsSize + MapSize + ApSize + AiSize + AxSize + 
                     14*dimDescSize + 6*ptrDescSize + 13*dimSuperSize + 3*ptrSuperSize +
                     2*nbatch*sizeof(int) + sizeof(int);


        /* case if exceed GPU memory */
        if(gpu_memtot >= Common->dev_mempool_size) {

  	/* store previous values */
          if(gpu_memtot_prev) {
            gpu_memtot = gpu_memtot_prev; 
            nbatch = nbatch_prev;
            maxsubtreeCsize = maxsubtreeCsize_prev;
            maxsubtreendesc = maxsubtreendesc_prev;
    	    maxsubtreebatch = maxsubtreebatch_prev;
            maxnumdescendantsperlevel = maxnumdescendantsperlevel_prev;
          }
          /* exit loop */
          break;
        }
        /* case if reached largest batch size in level */
        else if(nbatch == nsupernodes || nbatch >= MAXBATCHSIZE) {
          /* exit loop */
          break;
        }

    
        /* increment batch size */   
        nbatch_prev = nbatch;
        nbatch += 1;

      } /* end while loop */


      /* store max variables */
      if(maxsubtreeCsize > gb_p->maxCsize) gb_p->maxCsize = maxsubtreeCsize;		/* maximum C buffer size in any given subtree */
      if(maxsubtreendesc > gb_p->maxndesc) gb_p->maxndesc = maxsubtreendesc;		/* maximum # descendants in any given subtree */
      if(maxsubtreebatch > gb_p->maxbatch) gb_p->maxbatch = maxsubtreebatch;		/* maximum batch size in any given subtree */


    } 
    /* 
     * case if:
     * 1. CPU only  
     *
     */
    else if (runType == 1 || runType == 3) {    

      maxnumdescendantsperlevel = 0;
      nbatch = MAXBATCHSIZE;

      /* loop over supernodes in level */
      for(i = 0; i < nsupernodes; i++) {

        /* get supernode */                
        node = supernode_levels[supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+num_levels]+i];

        /* reset variables (new batch) */
        if( !(i % nbatch ) ) {
          Csize = 0;
          ndesc = 0;
        }

        Csize += supernode_size_desc[node];                 	/* add to total size of C buffer needed for storing schur complement in current batch of supernodes */
        ndesc += ndescendants[node];                        	/* add to total # descendants in current batch of supernodes */

        if(Csize > gb_p->maxCsize)  gb_p->maxCsize = Csize;           	/* store maximum C buffer size in any given level */
        if(ndesc > gb_p->maxndesc)  gb_p->maxndesc = ndesc;           	/* store maximum # descendants in any given level */
        if(nbatch > gb_p->maxbatch) gb_p->maxbatch = nbatch;          	/* store maximum batch size in any given level */
        if(ndesc > maxnumdescendantsperlevel) maxnumdescendantsperlevel = ndesc;
      } /* end loop over supernodes in level */
    }





    supernode_batch[supernode_levels_subtree_ptrs[subtree]+num_levels] = nbatch;	/* batch size per level */
    level_num_desc[count0++] = maxnumdescendantsperlevel;                     		/* total # descendants in current level */

    /* increment level */
    num_levels++;

    /* store pointer to level */
    supernode_levels_ptrs[supernode_levels_subtree_ptrs[subtree]+num_levels] = count2;

  } /* end loop over levels */





  /* store number of levels per subtree */
  supernode_num_levels[subtree] = num_levels;
      
  /* store counts */  
  counts[0] = count0;
  counts[1] = count1;
  counts[2] = count2;
  
}










/*
 *  Function:
 *    get_factor_size
 *
 *  Description:
 *    computes the size of the subfactor of the subtree 
 *
 */
void TEMPLATE2 (CHOLMOD (get_factor_size))
  (
   cholmod_global_pointers *gb_p,
   cholmod_cpu_pointers *cpu_p,
   cholmod_tree_pointers *tree_p,
   Int numSuper,
   Int subtree,
   Int *max_factor_size,
   Int *LpxSub
   )
{

  /* local variables */
  Int p, i, s, nscol, nsrow;
  Int *Super, *Lpi, *supernode_subtree, *supernode_subtree_ptrs, *factor_size;
  int numSubtree;

  /* set variables */
  p = 0;
  numSubtree	= gb_p->numSubtree;

  /* set host pointers */
  Super		= cpu_p->Super;
  Lpi		= cpu_p->Lpi;

  /* set tree pointers */
  supernode_subtree	= tree_p->supernode_subtree;
  supernode_subtree_ptrs	= tree_p->supernode_subtree_ptrs;
  factor_size		= tree_p->factor_size;



  /*
   * Store factor size in current subtree:  
   * For each subtree, calculate and store size of subfactor 
   * (Lxsub). Only do this for subtrees that go to GPU subtrees
   * algorithm (not last/root subtree). Also store largest
   * subfactor size, of all subtrees.
   */

  /* case:
   *   1. subtrees that go to GPU only algorithm (not last/root subtree)
   */
  if(subtree != numSubtree-1 /*|| numSubtree == 1*/)  {

    /* loop over supernodes */
    for(i=0; i < numSuper; i++) {

      /* get size of size of factor for each subtree */
      s = supernode_subtree[supernode_subtree_ptrs[subtree] + i];   
      nscol = Super [s+1] - Super [s] ;
      nsrow = Lpi[s+1] - Lpi[s] ;
      LpxSub [s] = p ;                                           /* store pointers to supernodes in sub-factor */
      p += nscol * nsrow ;                                       /* increment pointer to supernodes */
    } /* end loop over supernodes */
  } /* end case */



  /* store size of sub-factor for each subtree */
  factor_size[subtree] = p;

  /* store size of largest sub-factor of all subtrees */
  if(factor_size[subtree] > (*max_factor_size)) (*max_factor_size) = factor_size[subtree];

}










/*
 *  Function:
 *    gpu_num_descendants
 *
 *  Description:
 *    finds # descendants in supernode
 *
 */
void TEMPLATE2 (CHOLMOD (gpu_num_descendants))
  (
   cholmod_common *Common,
   cholmod_cpu_pointers *cpu_p,
   cholmod_tree_pointers *tree_p,
   Int s
   )
{

  Int d, nextd, n_descendant = 0;

  d = cpu_p->Head[s];
  while ( d != EMPTY )
  {
    nextd = cpu_p->Next[d];
    n_descendant++;
    d = nextd;
  }

  tree_p->ndescendants[s] = n_descendant;
}


/*
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
*/


