/*
 * File:
 *   cholmod_gpu
 *
 * Description:
 *   Contains functions for initializing 
 *   the GPU. 
 *
 */


/* includes */
#include "cholmod_internal.h"
#include "cholmod_core.h"
#include "cholmod_subtree.h"
#include "stdio.h"

#define MINSIZE (64 * 1024 * 1024)










/* 
 *  Function:
 *    poll_gpu
 *
 *  Description:
 *    Ensures a GPU is available. Returns true if it
 *    available, and false otherwise. 
 *
 */
static int poll_gpu (size_t s)          
{
#ifdef SUITESPARSE_CUDA
    /* local variables */
    void *p = NULL;


    /* Returns TRUE if the GPU has a block of memory of size s,
       FALSE otherwise. The block of memory is immediately freed. */
    if (s == 0) {
        return (FALSE) ;
    }
    if (cudaMalloc (&p, s) != cudaSuccess) {
        return (FALSE) ;
    }

    /* free block of memory */
    cudaFree (p) ;

    return (TRUE) ;
#endif
}










/*
 *  Function:
 *    cholmod_gpu_memorysize
 *
 *  Description:
 *     Determine the amount of free memory on the current GPU.  To use another
 *     GPU, use cudaSetDevice (k) prior to calling this routine, where k is an
 *     integer in the range 0 to the number of devices-1.   If the free size is
 *     less than 64 MB, then a size of 1 is returned. Function returns GPU 
 *     memory size available. 
 *     
 */
int CHOLMOD(gpu_memorysize)      
(
    size_t         *total_mem,
    size_t         *available_mem,
    cholmod_common *Common
)
{
    /* local variables */
    int k;
    size_t good, bad, s, total_free, total_memory;


    /* reset variables */
    *total_mem = 0;
    *available_mem = 0;


    /* early exit*/
#ifndef DLONG
    return 0;
#endif

    /* early exit */
    if (Common->useGPU != 1)
    {
        return (0) ;                    /* not using the GPU at all */
    }




    /* only if compiled with GPU */
#ifdef SUITESPARSE_CUDA


    /* find the total amount of free memory */
    cudaMemGetInfo (&total_free, &total_memory) ;
    *total_mem = total_memory;


    /* early exit (not even 64MB, return failure code)  */
    if (total_free < MINSIZE)
    {
        return (1) ;                    
    }


    /* try a bit less than the total free memory */
    s = MAX (MINSIZE, total_free*0.98) ;
    if (poll_gpu (s))
    {
        *available_mem = s;
        return (0) ;  
    }


    /* eary exit (not even 64MB, returnfailute code) */
    if (!poll_gpu (MINSIZE))
    {
        return (1) ;                    
    }


    /* 8 iterations of binary search */
    good = MINSIZE ;                    /* already known to be OK */
    bad  = total_free ;                 /* already known to be bad */
    for (k = 0 ; k < 8 ; k++)
    {
        s = (good + bad) / 2 ;
        if (poll_gpu (s))
        {
            good = s ;                  /* s is OK, increase good */
        }
        else
        {
            bad = s ;                   /* s failed, decrease bad */
        }
    }


    *available_mem = good;

#endif

    return (0) ; 
}










/*
 *  Function:
 *    cholmod_gpu_probe
 *
 *  Description:
 *   Used to ensure that a suitable GPU is available.  As this version of
 *   CHOLMOD can only utilize a single GPU, only the default (i.e. selected as
 *   'best' by the NVIDIA driver) is verified as suitable.  If this selection
 *   is incorrect, the user can select the proper GPU with the
 *   CUDA_VISIBLE_DEVICES environment variable.
 *
 *   To be considered suitable, the GPU must have a compute capability > 1 and
 *   more than 1 GB of device memory.
 */
int CHOLMOD(gpu_probe) ( cholmod_common *Common )
{
#ifdef SUITESPARSE_CUDA

    /* local variables */
    int ngpus, idevice, count=0, gpuid;
    double tstart, tend;
    struct cudaDeviceProp gpuProp, gpuProp0;




    /* early exit */
    if (Common->useGPU != 1)
    {
        return (0) ;
    }


    /* make sure requested # GPUs does not exceed max */
    if(Common->numGPU > CHOLMOD_MAX_NUM_GPUS)
    {
      Common->numGPU = CHOLMOD_MAX_NUM_GPUS;
    }


    /* make sure # GPUs does not exceed GPUs avaiable */
    cudaGetDeviceCount(&ngpus);
    if(Common->numGPU > ngpus)
    {
      Common->numGPU = ngpus;
    }
    


    /* if default selected (no numGPUs set) */
    if (Common->numGPU == -1) 
    {

      /* set GPU 0, fetch GPU 0, fetch GPU 0 properties */
      cudaSetDevice ( 0 );
      cudaGetDevice ( &idevice );
      cudaGetDeviceProperties ( &gpuProp0, idevice );
      count = 1;

      /* loop over all other GPUs */
      for(gpuid = 1; gpuid < ngpus; gpuid++) {

        /* set GPU, fetch GPU, fetch GPU properties  */
        cudaSetDevice ( gpuid );
        cudaGetDevice ( &idevice );
        cudaGetDeviceProperties ( &gpuProp, idevice );

        if(gpuProp.major == gpuProp0.major) {
          count++;
        }
                  
      }       

      /* set # GPUs */
      Common->numGPU = count;

    }
    /* if numGPUs set */
    else if (Common->numGPU > 0) 
    {

      /* loop over available gpus */
      for(gpuid = 0; gpuid < ngpus; gpuid++) 
      {
        
        /* set GPU, fetch GPU, fetch GPU properties  */
        cudaSetDevice ( gpuid );
        cudaGetDevice ( &idevice );
        cudaGetDeviceProperties ( &gpuProp, idevice );


        /* GPU must have at least:
         *   1. compute capability > 1
         *   2. > 1GB device memory */
        if ( gpuProp.major > 1 && 1.0e-9*gpuProp.totalGlobalMem > 1.0 )
            count++;				/* increment # qualified GPUs */        

      }

      /* reset # GPUs to # compatible devices */
      if(count < Common->numGPU)
        Common->numGPU = count;
            
    }


    if(Common->numGPU >= 1) return 1;		/* at least 1 compatible device found */
    else		    return 0;		/* no compatible devices, disable GPU */

#endif


    /* no GPU is available */
    return 0;  
}










/*
 *  Function:
 *    gpu_deallocate
 *
 *  Description:
 *    deallocate all GPU related buffers
 */
int CHOLMOD(gpu_deallocate)
(
    int gpuid,
    cholmod_common *Common
)
{
#ifdef SUITESPARSE_CUDA

    /* local variables */
    cudaError_t cudaErr;



    /* free device memory */
    if ( Common->dev_mempool[gpuid] )
    {
        cudaErr = cudaFree (Common->dev_mempool[gpuid]);
        if ( cudaErr ) {
            ERROR ( CHOLMOD_GPU_PROBLEM, "GPU error when freeing device memory.");
        }
    }


    Common->dev_mempool[gpuid] = NULL;
    Common->dev_mempool_size = 0;




    /* free host (pinned) memory */
    if ( Common->host_pinned_mempool[gpuid] )
    {
        cudaErr = cudaFreeHost ( Common->host_pinned_mempool[gpuid] );
        if ( cudaErr ) {
            ERROR ( CHOLMOD_GPU_PROBLEM, "GPU error when freeing host pinned memory.");
        }
    }


    Common->host_pinned_mempool[gpuid] = NULL;
    Common->host_pinned_mempool_size = 0;


  
    /* gpu end */
    CHOLMOD (gpu_end) (Common) ;
#endif

    return (0);
}










/* 
 *  Function:
 *    gpu_end
 *
 *  Description:
 *    free GPU handles & streams  
 *
 */
void CHOLMOD(gpu_end)
(
    cholmod_common *Common
)
{
#ifdef SUITESPARSE_CUDA

    /* local variables */
    int j, k;


    /* 
     * destroy cuBlas/cuSolver handles  
     */
    for(k = 0; k < Common->numGPU; k++) 
    {

      /* cuBlas handle */
      if (Common->cublasHandle[k])
      {
          cublasDestroy (Common->cublasHandle[k]) ;
          Common->cublasHandle[k] = NULL ;
      }


      /* cuSolver handle */
      if (Common->cusolverHandle[k])
      {
          cusolverDnDestroy (Common->cusolverHandle[k]) ;
          Common->cusolverHandle[k] = NULL ;
      }

    }




    /* 
     * destroy CUDA streams 
     */
    for(k = 0; k < Common->numGPU; k++) 
    {
      for(j = 0; j < CHOLMOD_DEVICE_STREAMS; j++) 
      {

        /* CUDA streams */
        if (Common->gpuStream[k][j])
        {
            cudaStreamDestroy (Common->gpuStream[k][j]) ;
            Common->gpuStream[k][j] = NULL ;        
        }

      }
    }




    /* 
     * destroy CUDA events 
     */
    for(k = 0; k < Common->numGPU; k++)
    {


      /* updateCKernelsComplete event */
      if (Common->updateCKernelsComplete[k])
      {
          cudaEventDestroy (Common->updateCKernelsComplete[k]) ;
          Common->updateCKernelsComplete[k] = NULL ;
      }


      /* updateCBuffersFree event */
      for(j = 0; j < CHOLMOD_HOST_SUPERNODE_BUFFERS; j++)
      {

        if (Common->updateCBuffersFree[k][j])
        {
            cudaEventDestroy (Common->updateCBuffersFree[k][j]) ;
            Common->updateCBuffersFree[k][j] = NULL ;
        }

      }


      /* cublasEventPotrf event */
      for(j = 0; j < 3; j++)
      {

        if (Common->cublasEventPotrf[k][j])
        {
            cudaEventDestroy (Common->cublasEventPotrf[k][j]) ;
            Common->cublasEventPotrf[k][j] = NULL ;
        }

      }
    }


#endif
}










/*
 *  Function: 
 *    gpu_allocate
 *
 *  Description:
 *    Allocate both host (pinned) and device memory needed for GPU computation.
 *    Cleans up the memory (memset, cudaMemset) to avoid linering errors/bad values
 *    from previous runs. Also creates cuBlas & cuSolver handles.
 *
 */

int CHOLMOD(gpu_allocate) 
( 
  cholmod_common *Common 
)
{

/* only if compiled with GPU */
#ifdef SUITESPARSE_CUDA


    /* local variables */
    int i, k;
    double tstart, tend;
    size_t fdm, tdm, tfdm[Common->numGPU], availableDeviceMemory, availableHostMemory, requestedDeviceMemory, requestedHostMemory, maxGpuMemBytes;
    cudaError_t cudaErr;
    cublasStatus_t cublasErr;
    cusolverStatus_t cusolverErr;


    /* early exit */
    if (Common->useGPU != 1) return (0) ;





    /* ensure valid input */
    maxGpuMemBytes = Common->maxGpuMemBytes;
    if ( maxGpuMemBytes < 0 ) maxGpuMemBytes = 0;






    /* 
     * Fetch total available device memory
     */
    fdm = 0;
    for(k = 0; k < Common->numGPU; k++) {

      cudaSetDevice(k); 
      CHOLMOD_HANDLE_CUDA_ERROR (CHOLMOD(gpu_memorysize) (&tdm,&tfdm[k],Common), "gpu_memorysize");

      /* get minimum amount of memory avialble across the 4 devices. Amount allocated for each 
       * device must be the same. */
      if(!k) fdm = tfdm[k];
      else {
        if(tfdm[k] < fdm) fdm = tfdm[k];
      }     

    }





    /* 
     * Compute the amount of device & host memory available: 
     * Always leave 50 MB free for driver use. 
     */
    availableDeviceMemory = fdm + Common->dev_mempool_size - 1024ll*1024ll*50ll;
    availableHostMemory = availableDeviceMemory;





    /* if memory requested larger than total memory available or no specific memory requested */
    if ( maxGpuMemBytes > availableDeviceMemory || maxGpuMemBytes == 0 ) {		
       maxGpuMemBytes = availableDeviceMemory;
    }





    /* 
     * Compute the amount of device & host memory requested
     */
    requestedDeviceMemory = maxGpuMemBytes;
    requestedHostMemory = maxGpuMemBytes;





    /* Set maximum size of buffers */
    Common->devBuffSize = requestedDeviceMemory/(size_t)(CHOLMOD_DEVICE_SUPERNODE_BUFFERS);
    Common->devBuffSize -= Common->devBuffSize%0x20000;





    /* deallocate memory for each GPU */
    for(k = 0; k < Common->numGPU; k++) {
      cudaSetDevice(k);
      CHOLMOD(gpu_deallocate) (k, Common);
    }





    /* allocated and clear pinned host memory for each GPU */
    for(k = 0; k < Common->numGPU; k++) {

      cudaSetDevice(k);

      /* allocate pinned memory */
      cudaErr = cudaHostAlloc ((void**)&(Common->host_pinned_mempool[k]), requestedHostMemory, cudaHostAllocMapped);    

      /* clear pinned memory */
      memset(Common->host_pinned_mempool[k],0,requestedHostMemory);

      if(cudaErr) printf("error:%s\n",cudaGetErrorString(cudaErr) );
      CHOLMOD_HANDLE_CUDA_ERROR (cudaErr,"host memory allocation failure\n");

    }





    /* allocate device memory for each GPU */
    for(k = 0; k < Common->numGPU; k++) {

      cudaSetDevice(k);

      /* allocate device memory */
      cudaErr = cudaMalloc ( &(Common->dev_mempool[k]), requestedDeviceMemory );

      /* clear device memory */
      cudaMemset(Common->dev_mempool[k],0,requestedDeviceMemory);

      if(cudaErr) printf("error:%s\n",cudaGetErrorString(cudaErr) );
      CHOLMOD_HANDLE_CUDA_ERROR (cudaErr,"device memory allocation failure\n");
    }





    /* store device & host memory sizes */
    Common->host_pinned_mempool_size = requestedHostMemory;
    Common->dev_mempool_size = requestedDeviceMemory;





    /* print GPU info */    
    PRINTF("\n\nGPU allocate..\n");
    PRINTFV("numGPU: %d\t",Common->numGPU);
    PRINTFV("maxGpuMemBytes: %ld\n",maxGpuMemBytes);
    PRINTFV("devBuffSize: %ld\n",Common->devBuffSize);
    PRINTFV("availableDeviceMemory:% ld\t",availableDeviceMemory);
    PRINTFV("availableHostMemory: %ld\n",availableHostMemory);
    PRINTFV("requestedDeviceMemory:% ld\t",requestedDeviceMemory);
    PRINTFV("requestedHostMemory: %ld\n",requestedHostMemory);
    PRINTF("\n\n\n");





    /* Create CUDA streams */
    for(k = 0; k < Common->numGPU; k++) {

      cudaSetDevice(k);

      for (i = 0; i < CHOLMOD_DEVICE_STREAMS; i++ ) {
        cudaErr = cudaStreamCreate ( &(Common->gpuStream[k][i]) );
        /* commenting this in causes occasinal nan values for benchmarks: 2cubes_sphere.mtx, crankseg_1.mtx */
        /*cudaErr = cudaStreamCreateWithFlags ( &(Common->gpuStream[k][i]), cudaStreamNonBlocking);*/
        if (cudaErr != cudaSuccess) {
          ERROR (CHOLMOD_GPU_PROBLEM, "create CUDA streams") ;
          abort();
        }
      }

    }





    /* Create CUDA handles */
    for(k = 0; k < Common->numGPU; k++) {

      cudaSetDevice(k);
      /* create cuBlas handle */
      cublasErr = cublasCreate (&(Common->cublasHandle[k])) ;
      if (cublasErr != CUBLAS_STATUS_SUCCESS) {
        ERROR (CHOLMOD_GPU_PROBLEM, "CUBLAS initialization") ;
        return 1;
      }

      /* create cuSolver handle */
      cusolverErr = cusolverDnCreate( &(Common->cusolverHandle[k]) );
      if (cusolverErr != CUSOLVER_STATUS_SUCCESS) {
        ERROR (CHOLMOD_GPU_PROBLEM, "CUSOLVER initialization") ;
        return 1;
      }

    }





    /* Create CUDA events */
    for(k = 0; k < Common->numGPU; k++) {

      cudaSetDevice(k);

      for (i = 0 ; i < 3 ; i++) {
        cudaErr = cudaEventCreateWithFlags(&(Common->cublasEventPotrf[k][i]), cudaEventDisableTiming) ;
        if (cudaErr != cudaSuccess) {
          ERROR (CHOLMOD_GPU_PROBLEM, "CUDA cublasEventPotrf event") ;
          return (0) ;
        }
      }

      for (i = 0 ; i < CHOLMOD_HOST_SUPERNODE_BUFFERS ; i++) {
        cudaErr = cudaEventCreateWithFlags(&(Common->updateCBuffersFree[k][i]), cudaEventDisableTiming) ;
        if (cudaErr != cudaSuccess) {
          ERROR (CHOLMOD_GPU_PROBLEM, "CUDA updateCBuffersFree event") ;
          return (0) ;
        }
      }

      cudaErr = cudaEventCreateWithFlags ( &(Common->updateCKernelsComplete[k]), cudaEventDisableTiming );
      if (cudaErr != cudaSuccess) {
        ERROR (CHOLMOD_GPU_PROBLEM, "CUDA updateCKernelsComplete event") ;
        return (0) ;
      }

    }





#endif
    return (0);
}
