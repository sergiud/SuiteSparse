/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_util.h ==================================== */
/* TODO check 2nd line of all files ... */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/*
 * General utility definitions for band reduction. None of these depend on
 * whether PIROBAND_FLOAT, PIROBAND_COMPLEX or PIROBAND_LONG is defined or not.
 */

/* This file should not be #include'd in user programs. */

#ifndef PIRO_BAND_UTIL_H
#define PIRO_BAND_UTIL_H

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>

/* ---------------------- Definition to turn on/off debugging -------------- */
/* Ensure that debugging is turned off */
#ifndef NDEBUG
#define NDEBUG
#endif

/* Uncomment/Comment next line to turn on/off debugging respectively */
/* #undef NDEBUG */

/* ---------------------- Definitions for MATLAB interfaces --------------- */
#ifdef MATLAB_MEX_FILE

/* TODO remove MALLOC and FREE */

#include "matrix.h"
#include "mex.h"
#define MALLOC mxMalloc
#define FREE mxFree

#ifndef NDEBUG
#define ASSERT(a) mxAssert(a, "")
#else
#define ASSERT(a)
#endif

#else /* MATLAB_MEX_FILE */

#ifndef NDEBUG
#include <assert.h>
#define ASSERT(a) assert(a)
#else
#define ASSERT(a)
#endif

#ifndef MALLOC
#define MALLOC malloc
#endif

#ifndef FREE
#define FREE free
#endif

#endif /* MATLAB_MEX_FILE */


/* Definition for the endcol() of a given row in the band matrix */
#define ENDCOL(k) ((k+bu > n-1) ? (n-1) : (k+bu))

#define SWAP(a, b, swap) { swap = (a) ; \
                        (a) = (b) ; \
                        (b) = swap ; \
                        } \

/* ========================================================================== */
/* ======== Utilities to ASSERT and PRINT in DEBUG mode ===================== */
/* ========================================================================== */

#ifndef NDEBUG  /* NDEBUG */

#define PRINT(params) printf params

#else /* NDEBUG */

#define PRINT(params)

#endif /* NDEBUG */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* Structure for safe minimum and safe maximum to compute the Givens
 * rotations reliably. They are computed once per matrix and used when
 * generating the Givens rotations.
 */

typedef struct
{
    double safemin2 ;
    double safemx2 ;
} piro_band_common ;

typedef struct
{
    float safemin2 ;
    float safemx2 ;
} piro_band_common_f ;

#endif /* PIRO_BAND_UTIL_H */
