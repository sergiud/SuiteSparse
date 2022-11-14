// =============================================================================
// === GPUQREngine/Include/GPUQREngine_Sentry.hpp ==============================
// =============================================================================

// GPUQREngine, Copyright (c) 2013, Timothy A Davis, Sencer Nuri Yeralan,
// and Sanjay Ranka.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------
//
// The SEntry struct is a tuple which is used to place the numeric values from
// the input problem into the dense representation of the fronts on the GPU.
//
// =============================================================================

#ifndef GPUQRENGINE_SENTRY_HPP
#define GPUQRENGINE_SENTRY_HPP

struct SEntry
{
    long findex;        // Index into the front where the value belongs
    double value;       // The numeric value to be placed

    SEntry()
    {
        findex = 0;
        value = 0.0;
    }
};

#endif
