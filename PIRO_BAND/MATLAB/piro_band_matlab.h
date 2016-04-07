/* ========================================================================== */
/* === PIRO_BAND/MATLAB/piro_band_matlab.h ================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Macro definitions for complex/real usage.
 * */

#ifndef PIRO_BAND_MATLAB_H
#define PIRO_BAND_MATLAB_H

/* maximum block size.  MAXBLK*MAXBLK must not cause integer overflow */
/* 64-by-64 is the typical maximum, so 2^15 is ample */
#define MAXBLK 32768

#undef MULT_I
#undef MULT_ADD_I
#undef MULT_SUB_I
#undef PRINT_QRVALUES

/* Make sure we use the PIROBAND_LONG version */
#define PIROBAND_LONG

#include "piro_band_util.h"
#include <string.h>

/* Include the appropriate MACROS for PIROBAND_LONG */
#include "piro_band_internal.h"
#include "piro_band_lapack_internal.h"


/* ========================================================================== */
/* ========================== Macros for Complex/Real usage ================= */
/* ========================================================================== */

#ifndef PIROBAND_COMPLEX
/* [ */

#define MULT_I(c, i1, a, i2, b, i3) (c[i1]) = (a[i2]) * (b[i3]) ;

#define MULT_ADD_I(c, i1, a, i2, b, i3) (c[i1]) += (a[i2]) * (b[i3]) ;

#define MULT_SUB_I(c, i1, a, i2, b, i3) (c[i1]) -= (a[i2]) * (b[i3]) ;

/* ] */
#else
/* [ */

#define MULT_I(c, i1, a, i2, b, i3)  \
                        (c[i1]) = (a[i2]) * (b[i3]) - (a[i2+1]) * (b[i3+1]) ; \
                        (c[i1+1]) = (a[i2]) * (b[i3+1]) + (a[i2+1]) * (b[i3]) ;

#define MULT_ADD_I(c, i1, a, i2, b, i3)  \
                        (c[i1]) += (a[i2]) * (b[i3]) - (a[i2+1]) * (b[i3+1]) ; \
                        (c[i1+1]) += (a[i2]) * (b[i3+1]) + (a[i2+1]) * (b[i3]) ;

#define MULT_SUB_I(c, i1, a, i2, b, i3)  \
                        (c[i1]) -= (a[i2]) * (b[i3]) - (a[i2+1]) * (b[i3+1]) ; \
                        (c[i1+1]) -= (a[i2]) * (b[i3+1]) + (a[i2+1]) * (b[i3]) ;

/* ] */
#endif

/* ========================================================================== */
/* ========================== Macros for printing complex/real data ========= */
/* ========================================================================== */
#ifndef NDEBUG

#ifdef PIROBAND_COMPLEX
#define PRINT_QRVALUES(str, i, val, i2) \
            mexPrintf(str "%0.4f %0.4f \n", i, val[i2+0], val[i2+1])
#else
#define PRINT_QRVALUES(str, i, val, i2) \
            mexPrintf(str "%0.4f \n", i, val[i2+0])
#endif

#else

#define PRINT_QRVALUES(str, i, val, i2)

#endif

/* ========================================================================== */
/* === options struct ======================================================= */
/* ========================================================================== */

#ifndef PIRO_BAND_MX_STRUCT
#define PIRO_BAND_MX_STRUCT

typedef struct piro_band_mx_options_struct
{
    Int blks [4] ;      /* blocksizes for block band reduction */

    Int sym ;           /* true if A is considered to be symmetric, in which
                           case only the upper triangular part of A is accessed.
                           The lower triangular part of A is assumed to be
                           the transpose of the upper triangular part (assuming
                           uplo is 'U' for the LAPACK interface).
                           False (0) otherwise. */

    Int econ ;          /* true if econ SVD is being computed, false otherwise.
                            only used by piro_band_svd. */

    Int benchmark ;     /* true if flops are to be computed, false otherwise */

    char uplo ;         /* 'U' or 'L' */

} piro_band_mx_options ;

#endif

/* ========================================================================== */

void piro_bandmex_storeband
(
    Int *Ap,
    Int *Ai,
    double *Ax,
    double *Axi,
    Int m,
    Int n,
    double *Cx,
    Int bu,
    Int bl,
    Int sym,
    char uplo
) ;

void piro_bandmex_find_full_bandwidth
(
    Int m,
    Int n,
    double *Ax,
    Int *bl,
    Int *bu
) ;

void piro_bandmex_storeband_withzeroes_full
(
    double *Ax,
    double *Axi,
    Int m,
    Int n,
    double *Cx,
    Int ub,
    Int lb
) ;

void piro_bandmex_band_conjugate_transpose
(
    Int m,
    Int n,
    Int bl,
    Int bu,
    double *A,
    double *AT,
    int iscomplex
) ;

mxArray *piro_bandmex_create_bidiagonal
(
    Int m,
    Int n,
    Int kind,
    double *b1,
    double *b2,
    Int offset
) ;

Int piro_bandmex_multiply
(
    Int a,
    Int b
) ;

Int piro_bandmex_add
(
    Int a,
    Int b
) ;

void piro_bandmex_error (Int err) ;

mxArray *piro_bandmex_put_dense
(
    Int m,
    Int n,
    double **Xhandle,
    int iscomplex,
    int transpose
) ;

void piro_bandmex_get_opts
(
    /* inputs */
    const mxArray *pargin [ ],  /* pargin from mexFunction */
    int nargin,                 /* nargin form mexFunction */
    int oparg,                  /* pargin [oparg] is the first option arg */

    /* output */
    piro_band_mx_options *opts
) ;

void piro_bandmex_blocksizes
(
    Int rc,
    Int cc,
    Int bl,
    Int bu,
    Int wantuv,
    piro_band_mx_options *opts
) ;

mxArray *piro_bandmex_identity
(
    Int n
) ;

mxArray *piro_bandmex_scalar
(
    double xr,      /* real part */
    double xi       /* imaginary part */
) ;

void piro_bandmex_find_bandwidth
(
    Int m,
    Int n,
    Int *Ap,
    Int *Ai,
    Int *bl,
    Int *bu
) ;

void piro_band_computeQ_drl
(
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int bl,                 /* lower bandwidth                         */
    double *V,              /* householder vectors                     */
    Int ldv,                /* leading dimension of V                  */
    double *Beta,           /* accumulated right rotations             */
    Int mq,                 /* #rows in Q                              */
    Int nq,                 /* #columns in Q                           */
    double *Q,              /* o/p Q of the QR factorization           */
    Int ldq,                /* leading dimension of Q                  */
    double *work,           /* workspace of size xxxxxxx               */
    double *X1              /* workspace of size xxxxxxx               */
) ;

void piro_band_computeQ_dcl
(
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int bl,                 /* lower bandwidth                         */
    double *V,              /* householder vectors                     */
    Int ldv,                /* leading dimension of V                  */
    double *Beta,           /* accumulated right rotations             */
    Int mq,                 /* #rows in Q                              */
    Int nq,                 /* #columns in Q                           */
    double *Q,              /* o/p Q of the QR factorization           */
    Int ldq,                /* leading dimension of Q                  */
    double *work,           /* workspace of size xxxxxxx               */
    double *X1              /* workspace of size xxxxxxx               */
) ;

int piro_band_qr_drl        /* returns error code */
(
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    double *A,              /* Band Matrix                             */
    Int ldab,               /* leading dimension of Ax                 */
    double *V,              /* o/p householder vectors                 */
    Int ldv,                /* leading dimension of V                  */
    double *Beta,           /* o/p accumulated right rotations         */
    double *work            /* workspace of size ((2 * bl) + bu + 1)   */
) ;

int piro_band_qr_dcl        /* returns error code */
(
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    double *A,              /* Band Matrix                             */
    Int ldab,               /* leading dimension of Ax                 */
    double *V,              /* o/p householder vectors                 */
    Int ldv,                /* leading dimension of V                  */
    double *Beta,           /* o/p accumulated right rotations         */
    double *work            /* workspace of size ((2 * bl) + bu + 1)   */
) ;

#endif
