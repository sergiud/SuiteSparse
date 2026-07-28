// =============================================================================
// === spqr_parallel ===========================================================
// =============================================================================

// SPQR, Copyright (c) 2008-2022, Timothy A Davis. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

// Factorize all the tasks in parallel with TBB.
// The GPU is not used.

#ifdef HAVE_TBB
#include "spqr.hpp"
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/task_group.h>

// =============================================================================
// === spqr_zippy ==============================================================
// =============================================================================

template <typename Entry, typename Int> SPQR_NO_EXPORT void spqr_zippy
(
    Int id,
    spqr_blob <Entry, Int> *Blob
)
{
    Int *TaskChildp = Blob->QRsym->TaskChildp ;
    Int *TaskChild  = Blob->QRsym->TaskChild ;
    Int pfirst = TaskChildp [id] ;
    Int plast  = TaskChildp [id+1] ;
    Int nchildren = plast - pfirst ;

    oneapi::tbb::task_group tasks ;
    for (Int i = 0 ; i < nchildren ; i++)
    {
        Int child = TaskChild [pfirst+i] ;
        tasks.run ([child, Blob] () {
            spqr_zippy <Entry, Int> (child, Blob) ;
        }) ;
    }
    tasks.wait () ;

    spqr_kernel <Entry, Int> (id, Blob) ;
}


// =============================================================================
// === spqr_parallel ===========================================================
// =============================================================================

template <typename Entry, typename Int> SPQR_NO_EXPORT void spqr_parallel
(
    Int ntasks,
    int nthreads,
    spqr_blob <Entry, Int> *Blob
)
{
    // Run the task tree, starting at the root id = ntasks-1.
    oneapi::tbb::task_arena arena (
        nthreads <= 0 ? oneapi::tbb::task_arena::automatic : nthreads) ;
    arena.execute ([ntasks, Blob] () {
        spqr_zippy <Entry, Int> (ntasks-1, Blob) ;
    }) ;
}
template void spqr_parallel <double, int32_t>
(
    int32_t ntasks,
    int nthreads,
    spqr_blob <double, int32_t> *Blob
) ;
template void spqr_parallel <Complex, int32_t>
(
    int32_t ntasks,
    int nthreads,
    spqr_blob <Complex, int32_t> *Blob
) ;
template void spqr_parallel <double, int64_t>
(
    int64_t ntasks,
    int nthreads,
    spqr_blob <double, int64_t> *Blob
) ;
template void spqr_parallel <Complex, int64_t>
(
    int64_t ntasks,
    int nthreads,
    spqr_blob <Complex, int64_t> *Blob
) ;
#endif
