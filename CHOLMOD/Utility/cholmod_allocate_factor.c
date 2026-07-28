//------------------------------------------------------------------------------
// CHOLMOD/Utility/cholmod_allocate_factor: allocate a simplicial factor
//------------------------------------------------------------------------------

// CHOLMOD/Utility Module. Copyright (C) 2023, Timothy A. Davis, All Rights
// Reserved.
// SPDX-License-Identifier: LGPL-2.1+

//------------------------------------------------------------------------------

// For backward compatibilty; L is returned as double precision.
// Use cholmod_alloc_factor to allocate a single precision factor.

#if defined (CHOLMOD_INT64)
#define CHOLMOD_ALLOCATE_FACTOR cholmod_l_allocate_factor
#define CHOLMOD_ALLOC_FACTOR cholmod_l_alloc_factor
#else
#if !defined (CHOLMOD_INT32)
#define CHOLMOD_INT32
#endif
#define CHOLMOD_ALLOCATE_FACTOR cholmod_allocate_factor
#define CHOLMOD_ALLOC_FACTOR cholmod_alloc_factor
#endif
#include "cholmod_internal.h"

cholmod_factor *CHOLMOD_ALLOCATE_FACTOR         // return the new factor L
(
    // input:
    size_t n,               // L is factorization of an n-by-n matrix
    cholmod_common *Common
)
{
    return (CHOLMOD_ALLOC_FACTOR (n, CHOLMOD_DOUBLE, Common)) ;
}
