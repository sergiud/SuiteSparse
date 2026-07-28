// =============================================================================
// === GPUQREngine/Source/LLBundle.cpp ================================
// =============================================================================

// GPUQREngine, Copyright (c) 2024, Timothy A Davis, Sencer Nuri Yeralan,
// and Sanjay Ranka.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------
//
// Instantiate LLBundle class template with types int32_t and int64_t
// for export from library.
//
// =============================================================================

#define GPUQRENGINE_NO_EXTERN_LLBUNDLE

#include "GPUQREngine_LLBundle.hpp"
#include "GPUQREngine_BucketList.hpp"

template class LLBundle<int32_t>;
template class LLBundle<int64_t>;
