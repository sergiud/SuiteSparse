/* ========================================================================== */
/* === Include/cholmod_subtree.h ============================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_subtree.h.
 * Copyright (C) 2016, Timothy A. Davis
 * CHOLMOD/Include/cholmod_subtree.h and the CHOLMOD GPU Module are licensed 
 * under Version 2.0 of the GNU General Public License.  See gpl.txt for a text 
 * of the license.  CHOLMOD is also available under other licenses; contact
 * authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */


#ifndef CHOLMOD_SUBTREE_H
#define CHOLMOD_SUBTREE_H

#ifdef SUITESPARSE_CUDA
#include "omp.h"
#include <fenv.h>
#endif




#define CHOLMOD_HANDLE_CUDA_ERROR(e,s) {if (e) {ERROR(CHOLMOD_GPU_PROBLEM,s);}}
#ifdef DLONG
#define Int SuiteSparse_long
#else
#define Int int
#endif




/* fine tuning parameters */
#define BINARY_SEARCH 5                         	/* # binary search steps to find ideal subtree size */
#define SUPERNODE_MIN 10                              	/* min. # supernodes in elimination tree to use subtrees alg. otherwise, everything sent to root tree */
#define MAXBATCHSIZE 100                               	/* max. # supernodes batched in a level (inversely proportional to subtree size) */
#define CHOLMOD_ND_ROW_LIMIT 256  			/* required descendant rows */
#define CHOLMOD_ND_COL_LIMIT 32   			/* required descendnat cols */
#define CHOLMOD_POTRF_LIMIT  512  			/* required cols for POTRF & TRSM on GPU */
#define CHOLMOD_GPU_SKIP     3  			/* # of host supernodes to perform before checking for free pinned buffers */ 




/* if verbose, enable prints & timers */
#ifdef CHOLMOD_VERBOSE
  /* GPU timers */
  #define TIMER_START(tstart,id)          	tstart[id]  = SuiteSparse_time()
  #define TIMER_END(tstart,tend,id)       	tend[id]   += SuiteSparse_time() - tstart[id]
  #define TIMER_START1(tstart)            	tstart      = SuiteSparse_time()
  #define TIMER_END1(tstart,tend,id)      	tend[id]   += SuiteSparse_time() - tstart
  #define CLEAR1(buf,id)                  	buf[id]     = 0
  #define SUM(buf,id,val)                 	buf[id]    += val
  /* prints */
  #define PRINTF(s)               		printf(s);
  #define PRINTFV(s,var)          		printf(s,var);
#else
  /* undefine GPU timers */
  #define TIMER_START(tstart,id)          
  #define TIMER_END(tstart,tend,id)       
  #define TIMER_START1(tstart)            
  #define TIMER_END1(tstart,tend,id)      
  #define CLEAR1(buf,id)                   
  #define SUM(buf,id,val)                 
  /* undefine prints */
  #define PRINTF(s)
  #define PRINTFV(s,var)
#endif




/* nvtx markers */
#ifdef SUITESPARSE_CUDA
#ifdef USE_NVTX
  #include "nvToolsExt.h"
#endif
#endif










/* struct for global variables */
typedef struct cholmod_global_pointers
{
    /* modes */
    int runType;    			/* determines type of run:	-1: CPU serial
					 *  				 0: GPU + CPU (hybrid)
					 *				 1: CPU parallel
					 *				 2: GPU only
					 *				 3: root only
					 */


    /* other global variables */
    int numGPU;				/* # GPUs to initialize ( < Common->useGPU if to few subtrees */
    int numDevice;                      /* # devices to use (GPU, CPU, root) */
    int numSubtree;			/* total # subtrees */
    int numRoot;			/* # root supernodes */
    int work_size;			/* size of workspace for cuSolver */
    int check[CHOLMOD_MAX_NUM_GPUS+2];	/* stores whether subtree is positive-defininte or not (1: positive-definite, 0: not positive definite) */
    
    /* max buffer sizes */
    Int maxCsize;			/* max Cbuff size (needed to store batch of schur complements) */
    Int maxndesc;			/* max # descendants in a batch */
    Int maxbatch;			/* max # supernodes in a batch */
    Int maxnsrow;			/* max nsrow (for cuSolver workspace) */
    Int maxnscol;			/* max nscol (for cuSolver workspace) */

    /* buffer sizes */
    size_t LxSize;                      /* size of factor 			(GPU) */
    size_t LsSize;			/* size of Ls 				(GPU) */
    size_t ApSize;			/* size of Ap 				(GPU) */
    size_t AiSize;			/* size of Ai				(GPU) */

    size_t AxSize;			/* size of Ax				(GPU) */
    size_t MapSize;			/* size of Map				(GPU & CPU) */
    size_t CSize;                       /* size of Cbuff			(GPU & CPU) */
    size_t dimDescSize;			/* size of dimensions array for desc.	(GPU) */
    size_t ptrDescSize;			/* size of pointers array for desc.	(GPU) */
    size_t dimSuperSize;		/* size of dimensions array for super	(GPU) */
    size_t ptrSuperSize;		/* size of pointers array for super	(GPU) */

    /* workspace buffers */
    Int *Iwork;					/* integer workspace */
    double *Xwork;				/* floating workspace */
    struct cholmod_subtree_order_t *Bwork;	/* subtree structure workspace */
    double *Cwork;                              /* workspace for Cbuff */
    Int *Mapwork;                               /* workspace for Map */

    /* workspace buffer size */
    size_t IworkSize;
    size_t XworkSize; 
    size_t BworkSize;	
    size_t CworkSize;
    size_t MapworkSize;

} cholmod_global_pointers ;










/* struct for load-balancing */
typedef struct cholmod_loadbalance_pointers
{
    /* load-balance arrays */
    Int    *numSubtreePerDevice;            /* # subtrees per device (GPU, CPU, root) */
    Int    *listSubtreePerDevice;           /* subtree id's per device  */
    double *workPerDevice;                 /* workload per device (measured in aggregate flop/flops) */
    double *subtreeSize;                    /* size of each subtree */
    struct cholmod_subtree_order_t *subtreeReorder;	 

} cholmod_loadbalance_pointers ;










/* struct for tree info */
typedef struct cholmod_tree_pointers
{
  /* tree arrays */
  Int *supernode_subtree;                        /* list of supernodes in each subtree      (pointer to supernode) */
  Int *supernode_subtree_ptrs;                   /* list of first supernode in each subtree (pointer to start of subtree) */
  Int *supernode_batch;
  Int *supernode_levels;                        /* list of supernodes in each level       (pointer to supernode) */
  Int *supernode_levels_ptrs;                   /* list of first supernode in each level  (pointer to start of level) */
  Int *supernode_levels_subtree_ptrs;            /* list of first supernode in each subtree (pointer to start of subtree) */
  Int *supernode_parent;                        /* list of parent in supernodes */
  Int *supernode_children;                      /* list of children in supernodes         (pointer to children) */
  Int *supernode_children_ptrs;                 /* list of first child in each supernode  (pointer to start of supernode) */
  Int *supernode_children_num;                  /* list of number of children in each supernode */
  Int *supernode_children_num2;                 /* list of number of children in each supernode - backup */
  Int *supernode_children_count;                
  Int *supernode_children_count2;               
  Int *supernode_num_levels;                    /* list of number of levels in each subtree */
  Int *level_descendants;                       /* list of maximum number of descendants any supernode has in each level */
  Int *level_descendants_ptrs;                  /* list of first (the above) in each subtree (pointer to start of subtree) */
  Int *level_num_desc;                          /* list of total number of descendants in each level */
  Int *level_num_desc_ptrs;                     /* list of first (the above) in each subtree (pointer to start of subtree) */
  Int *supernode_size_desc;                     /* list of total size of all descendants in each supernode */
  Int *supernode_size;                          /* list of size of each supernode */
  Int *ndescendants;                            /* list of number of descendants in each supernode */
  Int *supernode_root;                          /* list of supernode roots to each tree */
  Int *factor_size;                             /* list of sub-factor size for each subtree */
  double *supernode_flop;			/* list of floating-point operations (flop) per supernode for syrk,gemm,potrf,trsm operations */

} cholmod_tree_pointers ;










/* struct for CPU */
typedef struct cholmod_cpu_pointers
{

    /* host variables */
    Int Apacked;
    Int Fpacked;
    Int stype;
    double *beta;

    /* host buffers */
    Int *Ls; 
    Int *Lpi; 
    Int *Lpx; 
    Int *Lpos;
    Int *LpxSub;
    Int *Fp;
    Int *Fi;
    Int *Fnz;
    Int *Ap;
    Int *Ai;
    Int *Anz;
    Int *Super;
    Int *Map;
    Int *RelativeMap;
    Int *SuperMap;
    Int *Iwork;
    Int *Head;
    Int *Next;
    Int *Previous;
    Int *Next_save;
    Int *Lpos_save;
    double *Lx; 
    double *Ax;
    double *Az;
    double *Fx;
    double *Fz;
    double *C;

} cholmod_cpu_pointers ;










/* struct for GPU */
typedef struct cholmod_gpu_pointers
{
    int	   gpuid;

    /* device & pinned buffers */
    void   *gpuPtr[CHOLMOD_MAX_NUM_GPUS];
    void   *hostPtr[CHOLMOD_MAX_NUM_GPUS];
    double *h_Lx[CHOLMOD_MAX_NUM_GPUS];				/* factor in pinned memory */
    double *h_pLx[CHOLMOD_MAX_NUM_GPUS];			/* pointer to factor in pinned */
    double *d_Lx[CHOLMOD_MAX_NUM_GPUS];				/* factor in device */
    double *d_C[CHOLMOD_MAX_NUM_GPUS];				/* schur complement in device */
    Int    *d_Ls[CHOLMOD_MAX_NUM_GPUS];				/* Ls in device */
    Int    *d_Map[CHOLMOD_MAX_NUM_GPUS];			/* Map in device */    
    Int    *d_Ap[CHOLMOD_MAX_NUM_GPUS];				/* Ap in device */
    Int    *d_Ai[CHOLMOD_MAX_NUM_GPUS];				/* Ai in device */
    double *d_Ax[CHOLMOD_MAX_NUM_GPUS];				/* Ax in device */
    int    *d_info[CHOLMOD_MAX_NUM_GPUS];			/* info in device (potrf) */
    int    *d_devSync[CHOLMOD_MAX_NUM_GPUS]; 			/* devSync in device (potrf) */


    /* lists for batching supernodes */
    int    *h_dimSuper[CHOLMOD_MAX_NUM_GPUS];			/* ptr to supernode dim. in pinned */
    int    *d_dimSuper[CHOLMOD_MAX_NUM_GPUS];			/* ptr to supernode dim. in device */
    double **h_ptrSuper[CHOLMOD_MAX_NUM_GPUS];			/* ptr to supernode ptrs in pinned */
    double **d_ptrSuper[CHOLMOD_MAX_NUM_GPUS];			/* ptr to supernode ptrs in device */


    /* lists for batching descendants */
    int    *h_dimDesc[CHOLMOD_MAX_NUM_GPUS];   			/* ptr to descendant dim. in pinned */
    int    *d_dimDesc[CHOLMOD_MAX_NUM_GPUS];			/* ptr to descendant dim. in device */
    double **h_ptrDesc[CHOLMOD_MAX_NUM_GPUS];    		/* ptr to descendant ptrs in pinned */
    double **d_ptrDesc[CHOLMOD_MAX_NUM_GPUS];			/* ptr to descendant ptrs in device */


    /* structures for gemm,syrk,potrf,trsm parameters in device */
    struct cholmod_syrk_ptrs_t 	d_syrk[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_gemm_ptrs_t 	d_gemm[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_potrf_ptrs_t d_potrf[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_trsm_ptrs_t 	d_trsm[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_desc_ptrs_t 	d_desc[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_super_ptrs_t d_super[CHOLMOD_MAX_NUM_GPUS];

    /* structures for gemm,syrk,potrf,trsm parameters in host */
    struct cholmod_syrk_ptrs_t 	h_syrk[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_gemm_ptrs_t 	h_gemm[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_potrf_ptrs_t h_potrf[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_super_ptrs_t h_super[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_trsm_ptrs_t 	h_trsm[CHOLMOD_MAX_NUM_GPUS];
    struct cholmod_desc_ptrs_t 	h_desc[CHOLMOD_MAX_NUM_GPUS];

    
    /* root buffers */
    double *h_Lx_root[CHOLMOD_MAX_NUM_GPUS][CHOLMOD_HOST_SUPERNODE_BUFFERS] ;
    double *d_Lx_root[CHOLMOD_MAX_NUM_GPUS][2]; 
    double *d_C_root[CHOLMOD_MAX_NUM_GPUS];
    double *d_A_root[CHOLMOD_MAX_NUM_GPUS][2];
    Int *d_Ls_root[CHOLMOD_MAX_NUM_GPUS];
    Int *d_Map_root[CHOLMOD_MAX_NUM_GPUS];
    Int *d_RelativeMap_root[CHOLMOD_MAX_NUM_GPUS];

   
} cholmod_gpu_pointers ;










/* struct for profiling */
typedef struct cholmod_profile_pointers
{
    /* global timers */
    double g_start[20];
    double g_end[20];

    /* subtree timers */
    double b_start[CHOLMOD_MAX_NUM_GPUS+2];
    double b_end[CHOLMOD_MAX_NUM_GPUS+2];

    /* factorize timers */
    double f_start[CHOLMOD_MAX_NUM_GPUS+2][20];
    double f_end[CHOLMOD_MAX_NUM_GPUS+2][20];

    /* Blas timers */
    double syrk_time[CHOLMOD_MAX_NUM_GPUS+2][2];
    double gemm_time[CHOLMOD_MAX_NUM_GPUS+2][2];
    double potrf_time[CHOLMOD_MAX_NUM_GPUS+2][2];
    double trsm_time[CHOLMOD_MAX_NUM_GPUS+2][2];

    /* Blas flop */
    double syrk_flop[CHOLMOD_MAX_NUM_GPUS+2][2];
    double gemm_flop[CHOLMOD_MAX_NUM_GPUS+2][2];
    double potrf_flop[CHOLMOD_MAX_NUM_GPUS+2][2];
    double trsm_flop[CHOLMOD_MAX_NUM_GPUS+2][2]; 

} cholmod_profile_pointers ;










int cholmod_gpu_memorysize   /* GPU memory size available, 1 if no GPU */
(
    size_t         *total_mem,
    size_t         *available_mem,
    cholmod_common *Common
) ;

int cholmod_l_gpu_memorysize /* GPU memory size available, 1 if no GPU */
(
    size_t         *total_mem,
    size_t         *available_mem,
    cholmod_common *Common
) ;








 
int cholmod_gpu_probe   ( cholmod_common *Common ) ;
int cholmod_l_gpu_probe ( cholmod_common *Common ) ;

int cholmod_gpu_deallocate   ( int gpuid, cholmod_common *Common ) ;
int cholmod_l_gpu_deallocate ( int gpuid, cholmod_common *Common ) ;

void cholmod_gpu_end   ( cholmod_common *Common ) ;
void cholmod_l_gpu_end ( cholmod_common *Common ) ;

int cholmod_gpu_allocate   ( cholmod_common *Common ) ;
int cholmod_l_gpu_allocate ( cholmod_common *Common ) ;

#endif

