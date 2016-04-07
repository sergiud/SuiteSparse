/* ================== genband_util.h ==================================== */
/*
 * General utility definitions for band reduction. None of these depend on  
 * whether PIROBAND_FLOAT, PIROBAND_COMPLEX or PIROBAND_LONG is defined or not.
 */

#ifndef GENBAND_UTIL_H
#define GENBAND_UTIL_H

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>

/* ---------------------- Definition to turn on/off debugging -------------- */
/* Ensure that debugging is turned off */
#ifndef NDEBUG
#define NDEBUG
#endif

/* Ensure that printing is turned off */
#ifndef NPRINT
#define NPRINT
#endif

/* Uncomment/Comment next line to turn on/off debugging respectively */
/*
#undef NDEBUG
*/

/* Uncomment/Comment next line to turn on/off debug printing respectively */
/*
#undef NPRINT
*/

/* ---------------------- Definitions for MATLAB interfaces --------------- */

#ifdef MATLAB_MEX_FILE

#include "matrix.h"
#include "mex.h"
#define MALLOC mxMalloc
#define REALLOC mxRealloc
#define CALLOC mxCalloc
#define FREE mxFree

#ifndef NDEBUG
#define ASSERT(a, msg) mxAssert(a, msg)
#else
#define ASSERT(a, msg)
#endif

#else /* MATLAB_MEX_FILE */

#ifndef NDEBUG
#include <assert.h>
#define ASSERT(a, msg) assert(a)
#else
#define ASSERT(a, msg)
#endif

#ifndef MALLOC
#define MALLOC malloc
#endif

#ifndef FREE
#define FREE free
#endif

#ifndef REALLOC
#define REALLOC realloc
#endif

#ifndef CALLOC
#define CALLOC calloc
#endif

#endif /* MATLAB_MEX_FILE */

/* ---------------------- Internal Definition for benchmarking purposes ----- */
/* This definition is to count the number of floating point operations. (for 
 * internal use). If BENCHMARK is enabled when compiling the library then 
 * the global variable genband_flops has to be defined to use the library. This 
 * is the accurate flops for the real version of the library, not counting the
 * floating point operations for scaling while finding the givens rotations.
 * */
#ifdef BENCHMARK
extern double genband_flops ;

#define FLOPS(fl) genband_flops += ( ((fl) > 0) ? (fl): 0)

#else

#define FLOPS(fl)

#endif /* BENCHMARK */


/* Definition for the endcol() of a given row in the band matrix */
#define ENDCOL(k) ((k+bu > n-1) ? (n-1) : (k+bu))

#define SWAP(a, b, swap) { swap = (a) ; \
                        (a) = (b) ; \
                        (b) = swap ; \
                        } \

/* ========================================================================= */
/* ======= Utilities to ASSERT and PRINT in DEBUG mode ===================== */
/* ========================================================================= */

#ifndef NPRINT
/* enable diagnostic printing */
#define PRINT(params) printf params 
#define PRINT_INT_ARRAY(n,list,name) print_int_array (n, list, name)
#define PRINT_SKYLINE(sym,sky) print_skyline (sym, sky)
#else
/* no printing */
#define PRINT(params)
#define PRINT_INT_ARRAY(n,list,name)
#define PRINT_SKYLINE(sym,sky)
#endif

#if (defined(MATLAB_MEX_FILE) && !defined(NPRINT) && !defined(NDEBUG))
#define PRINT_SKYLINE_SVD(sym,sky,sgood) print_skyline_svd (sym,sky,sgood)
#else
#define PRINT_SKYLINE_SVD(sym,sky,sgood)
#endif

#ifndef NDEBUG  /* NDEBUG */

/* CHK_INDEX and CHK_CNT are there for debug purposes only. They are 
 * internal use. Do not uncomment */
/*#define CHK_INDEX(a, b, c, d) { ASSERT(INDEX((a), (b)) < 3100 && \
				INDEX((a), (b)) >= 0) ; \
				ASSERT(INDEX((c), (d)) < 3100 && \
				INDEX ((c), (d)) >= 0) ;  }

#define CHK_CNT(cnt) { ASSERT((cnt) >= 0 && (cnt) < 4) ; } */


#define CHK_CNT(cnt)
#define CHK_INDEX(a, b, c, d) 

#else /* NDEBUG */

#define CHK_CNT(cnt)
#define CHK_INDEX(a, b, c, d) 

#endif /* NDEBUG */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* ========================================================================= */
/* ======= Error codes for the entire package ======== ===================== */
/* ========================================================================= */
#define GENBAND_OK 0
#define GENBAND_BLKSIZE_INVALID (-1)
#define GENBAND_M_INVALID (-2)
#define GENBAND_N_INVALID (-3)
#define GENBAND_NRC_INVALID (-4)
#define GENBAND_BL_INVALID (-5)
#define GENBAND_BU_INVALID (-6)
#define GENBAND_AX_INVALID (-7)
#define GENBAND_LDAB_INVALID (-8)
#define GENBAND_B1_INVALID (-9)
#define GENBAND_B2_INVALID (-10)
#define GENBAND_U_INVALID (-11)
#define GENBAND_LDU_INVALID (-12)
#define GENBAND_V_INVALID (-13)
#define GENBAND_LDV_INVALID (-14)
#define GENBAND_C_INVALID (-15)
#define GENBAND_LDC_INVALID (-16)
#define GENBAND_WORK_INVALID (-17)
#define GENBAND_SYM_INVALID (-18)

/* ===========  Error codes for LAPACK style routines ======================= */

#define GENBAND_VECT_INVALID (-19)
#define GENBAND_OUT_OF_MEMORY (-20)
#define GENBAND_UPLO_INVALID (-21)


#endif /* GENBAND_UTIL_H */
