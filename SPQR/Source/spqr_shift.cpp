// =============================================================================
// === spqr_shift ==============================================================
// =============================================================================

// SPQR, Copyright (c) 2008-2022, Timothy A Davis. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

// Insert a zero as the first entry in a vector, shifting all the other entries
// down by one position.  Does nothing if X is NULL.

#include "spqr.hpp"

void spqr_shift
(
    // input, not modified
    int64_t n,

    // input/output
    int64_t *X                     // size n+1
)
{
    int64_t k ;
    if (X != NULL)
    {
        for (k = n ; k >= 1 ; k--)
        {
            X [k] = X [k-1] ;
        }
        X [0] = 0 ;
    }
}
