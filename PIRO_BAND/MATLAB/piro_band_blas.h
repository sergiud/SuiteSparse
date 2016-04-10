/* ========================================================================== */
/* === PIRO_BAND/MATLAB/piro_band_blas.h ==================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Macros for calling the correct LAPACK functions across various architectures.
 * */


#ifndef PIRO_BAND_BLAS_H
#define PIRO_BAND_BLAS_H

/* ========================================================================== */
/* === Architecture ========================================================= */
/* ========================================================================== */

#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define PIRO_BAND_SOL2
#define PIRO_BAND_ARCHITECTURE "Sun Solaris"

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define PIRO_BAND_SGI
#define PIRO_BAND_ARCHITECTURE "SGI Irix"

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define PIRO_BAND_LINUX
#define PIRO_BAND_ARCHITECTURE "Linux"

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define PIRO_BAND_AIX
#define PIRO_BAND_ARCHITECTURE "IBM AIX"
#define BLAS_NO_UNDERSCORE

#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define PIRO_BAND_ALPHA
#define PIRO_BAND_ARCHITECTURE "Compaq Alpha"

#elif defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#if defined (__MINGW32__) || defined (__MINGW32__)
#define PIRO_BAND_MINGW
#elif defined (__CYGWIN32__) || defined (__CYGWIN32__)
#define PIRO_BAND_CYGWIN
#else
#define PIRO_BAND_WINDOWS
#define BLAS_NO_UNDERSCORE
#endif
#define PIRO_BAND_ARCHITECTURE "Microsoft Windows"

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || defined (ARCH_HPUX)
#define PIRO_BAND_HP
#define PIRO_BAND_ARCHITECTURE "HP Unix"
#define BLAS_NO_UNDERSCORE

#elif defined (__hp700) || defined (MHP700) || defined (ARCH_HP700)
#define PIRO_BAND_HP
#define PIRO_BAND_ARCHITECTURE "HP 700 Unix"
#define BLAS_NO_UNDERSCORE

#else
/* If the architecture is unknown, and you call the BLAS, you may need to */
/* define BLAS_BY_VALUE, BLAS_NO_UNDERSCORE, and/or BLAS_CHAR_ARG yourself. */
#define PIRO_BAND_ARCHITECTURE "unknown"
#endif


/* ========================================================================== */
/* === BLAS and LAPACK names ================================================ */
/* ========================================================================== */

/* Prototypes for the various versions of the BLAS.  */

/* Determine if the 64-bit Sun Performance BLAS is to be used */
#if defined(PIRO_BAND_SOL2) && !defined(NSUNPERF) && defined(BLAS64)
#define SUN64
#endif

#ifdef SUN64

#define LAPACK_ZLARF zlarf_64_
#define LAPACK_DLARF dlarf_64_
#define LAPACK_ZLARFG zlarfg_64_
#define LAPACK_DLARFG dlarfg_64_
#define LAPACK_ZBDSQR zbdsqr_64_
#define LAPACK_DBDSQR dbdsqr_64_


#elif defined (BLAS_NO_UNDERSCORE)

#define LAPACK_ZLARF zlarf
#define LAPACK_DLARF dlarf
#define LAPACK_ZLARFG zlarfg
#define LAPACK_DLARFG dlarfg
#define LAPACK_ZBDSQR zbdsqr
#define LAPACK_DBDSQR dbdsqr

#else

#define LAPACK_ZLARF zlarf_
#define LAPACK_DLARF dlarf_
#define LAPACK_ZLARFG zlarfg_
#define LAPACK_DLARFG dlarfg_
#define LAPACK_ZBDSQR zbdsqr_
#define LAPACK_DBDSQR dbdsqr_

#endif


/* compile with -DBLAS64 if your LAPACK/BLAS uses 64-bit integers */
#if defined (LONGBLAS) || defined (BLAS64)
#define BLAS_INT SuiteSparse_long
#else
#define BLAS_INT int
#endif

/* TODO remove unused prototypes */

void LAPACK_DLARFT (char *direct, char *storev, BLAS_INT *n, BLAS_INT *k,
    double *V, BLAS_INT *ldv, double *Tau, double *T, BLAS_INT *ldt) ;

void LAPACK_ZLARFT (char *direct, char *storev, BLAS_INT *n, BLAS_INT *k,
    double *V, BLAS_INT *ldv, double *Tau, double *T, BLAS_INT *ldt) ;

void LAPACK_DLARFB (char *side, char *trans, char *direct, char *storev,
    BLAS_INT *m, BLAS_INT *n, BLAS_INT *k, double *V, BLAS_INT *ldv,
    double *T, BLAS_INT *ldt, double *C, BLAS_INT *ldc, double *Work,
    BLAS_INT *ldwork) ;

void LAPACK_ZLARFB (char *side, char *trans, char *direct, char *storev,
    BLAS_INT *m, BLAS_INT *n, BLAS_INT *k, double *V, BLAS_INT *ldv,
    double *T, BLAS_INT *ldt, double *C, BLAS_INT *ldc, double *Work,
    BLAS_INT *ldwork) ;

double BLAS_DNRM2 (BLAS_INT *n, double *X, BLAS_INT *incx) ;

double BLAS_DZNRM2 (BLAS_INT *n, double *X, BLAS_INT *incx) ;

void LAPACK_DLARFG (BLAS_INT *n, double *alpha, double *X, BLAS_INT *incx,
    double *tau) ;

void LAPACK_ZLARFG (BLAS_INT *n, double *alpha, double *X, BLAS_INT *incx,
    double *tau) ;

void LAPACK_DLARF (char *side, BLAS_INT *m, BLAS_INT *n, double *V,
    BLAS_INT *incv, double *tau, double *C, BLAS_INT *ldc, double *Work) ;

void LAPACK_ZLARF (char *side, BLAS_INT *m, BLAS_INT *n, double *V,
    BLAS_INT *incv, double *tau, double *C, BLAS_INT *ldc, double *Work) ;

void LAPACK_DBDSQR (char *uplo,
    BLAS_INT *n, BLAS_INT *ncvt, BLAS_INT *nru,  BLAS_INT *ncc,
    double *d, double *e, double *Vt, BLAS_INT *ldvt, double *U, BLAS_INT *ldu,
    double *C, BLAS_INT *ldc, double *work, BLAS_INT *info) ;

void LAPACK_ZBDSQR (char *uplo,
    BLAS_INT *n, BLAS_INT *ncvt, BLAS_INT *nru,  BLAS_INT *ncc,
    double *d, double *e, double *Vt, BLAS_INT *ldvt, double *U, BLAS_INT *ldu,
    double *C, BLAS_INT *ldc, double *work, BLAS_INT *info) ;

#endif
