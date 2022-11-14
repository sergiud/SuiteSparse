// =============================================================================
// === spqrgpu_computeFrontStaging =============================================
// =============================================================================

// SPQRGPU, Copyright (c) 2008-2022, Timothy A Davis, Sanjay Ranka,
// Sencer Nuri Yeralan, and Wissam Sid-Lakhdar, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

// Returns a front staging and whether the staging is feasible.
// A front staging is infeasible if a front and its children do not fit on
// the GPU at the same time.

#ifdef SUITESPARSE_CUDA
#include "spqr.hpp"
#include "GPUQREngine_Scheduler.hpp"

void spqrgpu_computeFrontStaging
(
    // inputs, not modified on output
    int64_t numFronts,     // total number of fronts (nf in caller)
    int64_t *Parent,       // size nf+1, assembly tree (f=nf is placeholder)
    int64_t *Childp,       // size nf+2, children of f are
                        //      Child [Childp [f] ... Childp [f+1]-1]
    int64_t *Child,        // size nf+1.

    int64_t *Fm,           // size nf+1, front f has Fm [f] rows
    int64_t *Cm,           // size nf+1, front f has Cm [f] rows in contrib
    int64_t *Rp,           // size nf+1, Rj[Rp[f]...Rp[f+1]-1] are the cols in f
    int64_t *Sp,           // size m+1, row pointers for sparse matrix S
    int64_t *Sleft,        // size n+2, see spqr_stranspose for description
    int64_t *Super,        // size nf+1, front f pivotal cols are
                        //      Super[f]..Super[f+1]-1
    int64_t *Post,         // size nf+1, front f is kth, if f = Post [k]

    int64_t RimapSize,     // scalar, size of Rimap on the GPU (# of int's)
    int64_t RjmapSize,     // scalar, size of Rimap on the GPU (# of int's)

    // output, not defined on input:
    bool *feasible,     // scalar, true if feasible, false if GPU memory too low
    int64_t *numStages,    // scalar, number of stages
    int64_t *Stagingp,     // size nf+2, fronts are in the list
                        //      Post [Stagingp [stage]...Stagingp[stage+1]-1]
    int64_t *StageMap,     // size nf, front f is in stage StageMap [f]

    size_t *FSize,      // size nf+1, FSize[stage]: size in bytes of MongoF
    size_t *RSize,      // size nf+1, Rsize[stage]: size in bytes of MongoR
    size_t *SSize,      // size nf+1, Ssize[stage]: size in bytes of S
    int64_t *FOffsets,     // size nf, front f in MondoF [FOffsets[f]...] on GPU
    int64_t *ROffsets,     // size nf, R block in MondoR [Roffsets[f]...] on CPU
    int64_t *SOffsets,     // size nf, S entries for front f are in
                        //      wsS [SOffsets[f]...]

    // input/output:
    cholmod_common *cc
)
{

    // -------------------------------------------------------------------------
    // determine available GPU memory, required for all stages
    // -------------------------------------------------------------------------

    // gpuMemorySize = 0 is for testing only.  This value forces each front to
    // appear in its own stage.  gpuMemorySize = 1 is also for testing, to
    // check how this function handles problem that is infeasible due to lack
    // of GPU memory.  We must also ensure that gpuMemorySize does not
    // accidentally become negative.

    size_t gpuMemorySize = cc->gpuMemorySize;

//  printf ("GPU mem starts %g MB\n", (double) gpuMemorySize / (1024*1024)) ;

    // account for two Scheduler work queues in the GPU memory
    if (gpuMemorySize > 1)
    {
        // The GPU must hold two workspace queues, each of size maxQueueSize,
        // and where each entry is of size sizeof(TaskDescriptor)
        size_t maxQueueSize = ssgpu_maxQueueSize (gpuMemorySize) ;
        size_t s = 2 * maxQueueSize * sizeof (TaskDescriptor) ;
        gpuMemorySize = (gpuMemorySize > s) ? (gpuMemorySize-s) : 0 ;
    }

    // account for Rimap in the GPU memory
    if (gpuMemorySize > 1)
    {
        size_t s = RimapSize * sizeof (int) ;
        gpuMemorySize = (gpuMemorySize > s) ? (gpuMemorySize-s) : 0 ;
    }

    // account for Rjmap in the GPU memory
    if (gpuMemorySize > 1)
    {
        size_t s = RjmapSize * sizeof (int) ;
        gpuMemorySize = (gpuMemorySize > s) ? (gpuMemorySize-s) : 0 ;
    }

    // account for cudaMalloc memory manager overhead in the GPU memory
    if (gpuMemorySize > 1)
    {
        size_t s = 1024 * 1024 ;        // just 1 MB for good measure
        gpuMemorySize = (gpuMemorySize > s) ? (gpuMemorySize-s) : 0 ;
    }

//  printf ("GPU mem now    %g MB\n", (double) gpuMemorySize / (1024*1024)) ;

    // -------------------------------------------------------------------------
    // assign fronts to stages based on remaining GPU memory availability
    // -------------------------------------------------------------------------

    /* The memory requirement for a front is the summation of the memory
       requirements from its children plus its own memory requirement.
       If we use a postorder traversal, we only need to keep a rolling sum. */
    size_t ReqMem = 0;   // also used in FOffsets

    /* RMem is always < ReqMem.
       RMem is not just for R, but an upper-bound on the amount of front data
       we have to pull back from the GPU in the event that the front is staged.
       We already would have room to pull back the CBlock if we account for it.
       The QREngine is intelligent enough to pull only the needed FrontData.
       The language is clearer in QREngine as well ("R" renamed "FrontData") */
    size_t RMem = 0;     // also used in ROffsets

    /* SMem is twice the number of S values (one for index, one for value). */
    size_t SMem = 0;     // also used in SOffsets

    /* VTMem is the amount of memory required for the VT blocks. */
    size_t VTMem = 0;

    int64_t stage = 0;
    Stagingp[0] = 0;
    for(int64_t p=0; p<numFronts; p++)
    {
        int64_t f = Post[p]; // The postordering ensures we visit children first

        int64_t fm = Fm[f];
        int64_t fn = Rp[f+1] - Rp[f];
        int64_t fp = Super[f+1] - Super[f];
        int64_t frank = MIN(fm, fp);
        // int64_t cn = fn - fp ;
        int64_t cm = Cm[f] ;
        size_t frontMem = fm * fn;          // F
        size_t rMem = (frank + cm) * fn;    // R + C

        // for sMem, "2 *" assumes sizeof (SEntry) is 2*sizeof(double)
        size_t sMem = 2 * (Sp[Sleft[Super[f+1]]] - Sp[Sleft[Super[f]]]);

        // CEIL is defined in GPUQREngine
        size_t vtMem = CEIL(fm, TILESIZE) * (TILESIZE+1) * TILESIZE ;

        size_t childMemOld = 0;
        size_t childMemNew = 0;

        /* Add contribution block memory from children in earlier stages. */
        for(int64_t cp=Childp[f]; cp<Childp[f+1]; cp++)
        {
            int64_t c = Child[cp];

            // int64_t cfm = Fm[c];
            int64_t cfn = Rp[c+1] - Rp[c];
            int64_t cfp = Super[c+1] - Super[c];
            // int64_t crank = MIN(cfm, cfp);
            int64_t ccn = cfn - cfp ;
            int64_t ccm = Cm[c];

            if(StageMap[c] < stage)
            {
                childMemOld += ccm * ccn;
            }
            else
            {
                childMemNew += ccm * ccn;
            }
        }

        /* determine which stage will contain this front */
        if((ReqMem + (frontMem + childMemOld + sMem + vtMem)) * sizeof(double)
            < gpuMemorySize)
        {
            /* If we can add the front to the current stage, accum its mem. */
            FOffsets[f] = ReqMem;
            ROffsets[f] = RMem;
            SOffsets[f] = SMem / 2; // correct for data width
            ReqMem += frontMem + childMemOld;
            RMem += rMem;
            SMem += sMem;
            VTMem += vtMem;
        }
        else if (gpuMemorySize == 0 ||
            ((frontMem + childMemOld + childMemNew + sMem + vtMem)
             * sizeof(double) < gpuMemorySize))
        {
            /* Else if the front and its children fit on the GPU, add it
               to the next stage and reset the mem accumulator. */
            PR (("gpuMemorySize: move front to next stage\n")) ;
            FSize[stage] = ReqMem;  // save the sizes for the stage
            RSize[stage] = RMem;
            SSize[stage] = SMem / 2; // correct for data width
            Stagingp[++stage] = p;
            ReqMem = 0;             // move onto the next stage
            RMem = 0;
            SMem = 0;
            VTMem = 0;
            FOffsets[f] = ReqMem;
            ROffsets[f] = RMem;
            SOffsets[f] = SMem / 2; // correct for data width
            ReqMem += frontMem + childMemOld + childMemNew;
            RMem += rMem;
            SMem += sMem;
            VTMem += vtMem;
        }
        else
        {
            /* Else the front and its children can't fit on the GPU,
               so we have an infeasible schedule. */
            PR (("gpuMemorySize too small: schedule infeasible\n")) ;
            ERROR (CHOLMOD_GPU_PROBLEM, "GPU memory too small\n") ;
            *numStages = 0 ;
            *feasible = false ;
            return ;
        }
        StageMap[f] = stage;
    }

    /* Make sure that even if everything fits in one stage
       that we finalize the stage. */
    FSize[stage] = ReqMem;
    RSize[stage] = RMem;
    SSize[stage] = SMem / 2;  // correct for data width
    Stagingp[++stage] = numFronts;
    if(Stagingp[stage] == Stagingp[stage-1]) stage--;

    *numStages = stage;
    *feasible = true;
    return;
}
#endif
