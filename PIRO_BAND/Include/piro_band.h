/* ========================================================================== */
/* === PIRO_BAND/Include/piro_band.h ======================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* External interface for the band reduction routines. There are eight different
 * versions for the reduce routine which reduces a band matrix A stored in
 * packed band format to a bidiagonal matrix such that A = U * B * V'
 * using Givens rotations. The left and the right transformations are
 * accumulated in U and V if they are present. Only the upper triangular part
 * of the symmetric matrices are expected when the sym flag is turned on.
 *
 * The functions in PIRO_BAND come in 8 different flavors, with all 8
 * combinations of double/single precision, real/complex, and int/"long".
 * A "long" is actually SuiteSparse_long, defined by SuiteSparse_config.h in
 * the SuiteSparse_config package.  It is "long" on most systems, except for
 * Windows 64 (where it is int64_t).
 *
 * The naming convention for the routines will be piro_band_reduce_<xyz>
 *  where x - d/s based on double/single precision.
 *        y - r/c based in real/complex matrices.
 *        z - i/p based on int/long.
 *  For example, the interface name for double, real and integer version will be
 *  piro_band_reduce_dri.
 * */

#ifndef PIRO_BAND_H
#define PIRO_BAND_H

/* make it easy for C++ programs to include piro_band_reduce */
#ifdef __cplusplus
extern "C" {
#endif

#include "SuiteSparse_config.h"

/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* double, real, int */
int piro_band_reduce_dri    /* returns PIRO_BAND_OK (0) if successful, a
                               negative error code otherwise (see error
                               codes below */
(
    int blks[],
    /* An array of size four for the block size in upper and lower bands.
     * Specifically :
     * blk[0], blk[1] for #columns, #rows in the block for upper band
     * blk[2], blk[3] for #rows, #columns in the block for lower band
     * Block sizes cannot be negative. They can be zero only if the
     * corresponding bandwidth is zero.
     * blk[0] + blk[1] <= upper bandwidth
     * blk[2] + blk[3] <= lower bandwidth + 1 if lower bandwidth > 0
     * blk[2] + blk[3] = 0 otherwise.
     * */
    int m,                  /* #rows in the original matrix.  m >= 0        */

    int n,                  /* #columns in the original matrix. n >= 0      */

    int nrc,                /* #rows in C
                             * nrc >= 0, and nrc = 0 if C == NULL           */

    int bl,                 /* lower bandwidth,
                             * bl >= 0 and bl <= m-1 if m > 0.              */

    int bu,                 /* upper bandwidth,
                             * bu >= 1 and bu <= n-1 if n > 0.              */

    double *A,              /* Band Matrix stored in packed band format,
    * If the matrix is complex then the real and imaginary parts should be
    * stored next to each other like C99 complex data type. zomplex where we
    * store real and imaginary parts are stored in two separate arrays are
    * not supported.
    * */

    int ldab,               /* leading dimension of A
                             * ldab >= (bl + bu + 1)
                             * If A is symmetric then ldab >= bu + 1        */

    double *B1,             /* Output : diagonal 0 of the bidiagonal matrix
                             * A vector of size min(m, n). Always real      */

    double *B2,             /* Output : diagonal 1 of the bidiagonal matrix
                             * A vector of size min(m, n)-1. Always real.   */

    double *U,              /* Output : accumalated left rotations, mxm
                             * Orthogonal matrix. The rotations are always
                             * accumulated, and U always initialized to
                             * identity if U != NULL                        */

    int ldu,                /* leading dimension of U
                             * ldu >= 0, especially ldu = 0 if U == NULL
                             * and ldu >= max(1, m) if U != NULL            */

    double *V,              /* Output :  accumalated right rotations, nxn
                             * orthogonal matrix. The rotations are always
                             * accumulated and V always initialized
                             * to identity if V != NULL.                    */

    int ldv,                /* leading dimension of V
                             * ldv >= 0, especially ldv = 0 if V == NULL
                             * and ldv >= max(1, n) if V != NULL            */

    double *C,              /* Output : for left rotations  nrc x m matrix
                             * which accumulates the left rotations. C is
                             * never initialized.
                             * */

    int ldc,                /* leading dimension of C, nrc >= 0, ldc >= nrc */

    double *dws,            /* workspace, to store the blocked Givens
    * rotations, of size 2*MAX(blk[0]*blk[1], blk[2]*blk[3]). If the input is
    * complex then we need twice this space to store the complex rotations.
    * In the complex case, We store both cosine and sine of the rotations as
    * complex scalars even though the cosine is always real for simpler
    * indexing.                                                             */

    int sym                 /* Flag : 0/1 to specify symmetric inputs       */
) ;

/* See comments in piro_band_reduce_dri above for detailed description of every
 * argument to this function                                                */
/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* double, complex, int */
int piro_band_reduce_dci    /* returns 0 if successful, < 0 on error */
(
    int blks[],
    int m,          /* #rows in the input matrix */
    int n,          /* #columns in the input matrix */
    int nrc,        /* #rows in C      */
    int bl,         /* lower bandwidth */
    int bu,         /* upper bandwidth */
    double *Ax,     /* input matrix in packed band format */
    int ldab,       /* leading dimension of A */
    double *B1,     /* output : the main diagonal of the bidiagonal matrix */
    double *B2,     /* output : the super diagonal of the bidiagonal matrix */
    double *U,      /* Accumulated left rotations */
    int ldu,        /* leading dimension of U */
    double *V,      /* Accumulated right rotations */
    int ldv,        /* leading dimension of V */
    double *C,      /*  Accumulated left rotations */
    int ldc,        /* leading dimension of C */
    double *dws,    /* workspace */
    int sym         /* 0/1 to indicate whether A is a symmetric matrix */
) ;


/* See comments in piro_band_reduce_dri above for detailed description of every
 * argument to this function                                                */
/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* single, real, int */
int piro_band_reduce_sri    /* returns 0 if successful, < 0 on error */
(
    int blks[],
    int m,          /* #rows in the input matrix */
    int n,          /* #columns in the input matrix */
    int nrc,        /* #rows in C      */
    int bl,         /* lower bandwidth */
    int bu,         /* upper bandwidth */
    float *Ax,      /* input matrix in packed band format */
    int ldab,       /* leading dimension of A */
    float *B1,      /* output : the main diagonal of the bidiagonal matrix */
    float *B2,      /* output : the super diagonal of the bidiagonal matrix */
    float *U,       /* Accumulated left rotations */
    int ldu,        /* leading dimension of U */
    float *V,       /* Accumulated right rotations */
    int ldv,        /* leading dimension of V */
    float *C,       /* Accumulated left rotations */
    int ldc,        /* leading dimension of C */
    float *dws,     /* workspace */
    int sym         /* 0/1 to indicate whether A is a symmetric matrix */
) ;

/* See comments in piro_band_reduce_dri above for detailed description of every
 * argument to this function                                                */
/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* single, complex, int */
int piro_band_reduce_sci    /* returns 0 if successful, < 0 on error */
(
    int blks[],
    int m,          /* #rows in the input matrix */
    int n,          /* #columns in the input matrix */
    int nrc,        /* #rows in C      */
    int bl,         /* lower bandwidth */
    int bu,         /* upper bandwidth */
    float *Ax,      /* input matrix in packed band format */
    int ldab,       /* leading dimension of A */
    float *B1,      /* output : the main diagonal of the bidiagonal matrix */
    float *B2,      /* output : the super diagonal of the bidiagonal matrix */
    float *U,       /* Accumulated left rotations */
    int ldu,        /* leading dimension of U */
    float *V,       /* Accumulated right rotations */
    int ldv,        /* leading dimension of V */
    float *C,       /* Accumulated left rotations */
    int ldc,        /* leading dimension of C */
    float *dws,     /* workspace */
    int sym         /* 0/1 to indicate whether A is a symmetric matrix */
) ;


/* See comments in piro_band_reduce_dri above for detailed description of every
 * argument to this function                                                */
/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* double, real, long */
int piro_band_reduce_drl    /* returns 0 if successful, < 0 on error */
(
    SuiteSparse_long blks[],
    SuiteSparse_long m,     /* #rows in the input matrix */
    SuiteSparse_long n,     /* #columns in the input matrix */
    SuiteSparse_long nrc,   /* #rows in C      */
    SuiteSparse_long bl,    /* lower bandwidth */
    SuiteSparse_long bu,    /* upper bandwidth */
    double *Ax,             /* input matrix in packed band format */
    SuiteSparse_long ldab,  /* leading dimension of A */
    double *B1,             /* output : main diagonal of bidiagonal matrix */
    double *B2,             /* output : super diagonal of bidiagonal matrix */
    double *U,              /* Accumulated left rotations */
    SuiteSparse_long ldu,   /* leading dimension of U */
    double *V,              /* Accumulated right rotations */
    SuiteSparse_long ldv,   /* leading dimension of V */
    double *C,              /* Accumulated left rotations */
    SuiteSparse_long ldc,   /* leading dimension of C */
    double *dws,            /* workspace */
    SuiteSparse_long sym    /* 0/1 indicates whether A is a symmetric matrix */
) ;

/* See comments in piro_band_reduce_dri above for detailed description of every
 * argument to this function                                                */
/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* double, complex, long */
int piro_band_reduce_dcl    /* returns 0 if successful, < 0 on error */
(
    SuiteSparse_long blks[],
    SuiteSparse_long m,     /* #rows in the input matrix */
    SuiteSparse_long n,     /* #columns in the input matrix */
    SuiteSparse_long nrc,   /* #rows in C      */
    SuiteSparse_long bl,    /* lower bandwidth */
    SuiteSparse_long bu,    /* upper bandwidth */
    double *Ax,             /* input matrix in packed band format */
    SuiteSparse_long ldab,  /* leading dimension of A */
    double *B1,             /* output : main diagonal of bidiagonal matrix */
    double *B2,             /* output : super diagonal of bidiagonal matrix */
    double *U,              /* Accumulated left rotations */
    SuiteSparse_long ldu,   /* leading dimension of U */
    double *V,              /* Accumulated right rotations */
    SuiteSparse_long ldv,   /* leading dimension of V */
    double *C,              /* Accumulated left rotations */
    SuiteSparse_long ldc,   /* leading dimension of C */
    double *dws,            /* workspace */
    SuiteSparse_long sym    /* 0/1 indicates whether A is a symmetric matrix */
) ;


/* See comments in piro_band_reduce_dri above for detailed description of every
 * argument to this function                                                */
/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* single, real, long */
int piro_band_reduce_srl    /* returns 0 if successful, < 0 on error */
(
    SuiteSparse_long blks[],
    SuiteSparse_long m,     /* #rows in the input matrix */
    SuiteSparse_long n,     /* #columns in the input matrix */
    SuiteSparse_long nrc,   /* #rows in C      */
    SuiteSparse_long bl,    /* lower bandwidth */
    SuiteSparse_long bu,    /* upper bandwidth */
    float *Ax,              /* input matrix in packed band format */
    SuiteSparse_long ldab,  /* leading dimension of A */
    float *B1,              /* output : main diagonal of bidiagonal matrix */
    float *B2,              /* output : super diagonal of bidiagonal matrix */
    float *U,               /* Accumulated left rotations */
    SuiteSparse_long ldu,   /* leading dimension of U */
    float *V,               /* Accumulated right rotations */
    SuiteSparse_long ldv,   /* leading dimension of V */
    float *C,               /* Accumulated left rotations */
    SuiteSparse_long ldc,   /* leading dimension of C */
    float *dws,             /* workspace */
    SuiteSparse_long sym    /* 0/1 indicates whether A is a symmetric matrix */
) ;

/* See comments in piro_band_reduce_dri above for detailed description of every
 * argument to this function                                                */
/* Prototypes for bidiagonal reduction of band matrices with the datatypes : */
/* single, complex, long */
int piro_band_reduce_scl    /* returns 0 if successful, < 0 on error */
(
    SuiteSparse_long blks[],
    SuiteSparse_long m,     /* #rows in the input matrix */
    SuiteSparse_long n,     /* #columns in the input matrix */
    SuiteSparse_long nrc,   /* #rows in C      */
    SuiteSparse_long bl,    /* lower bandwidth */
    SuiteSparse_long bu,    /* upper bandwidth */
    float *Ax,              /* input matrix in packed band format */
    SuiteSparse_long ldab,  /* leading dimension of A */
    float *B1,              /* output : main diagonal of bidiagonal matrix */
    float *B2,              /* output : super diagonal of bidiagonal matrix */
    float *U,               /* Accumulated left rotations */
    SuiteSparse_long ldu,   /* leading dimension of U */
    float *V,               /* Accumulated right rotations */
    SuiteSparse_long ldv,   /* leading dimension of V */
    float *C,               /* Accumulated left rotations */
    SuiteSparse_long ldc,   /* leading dimension of C */
    float *dws,             /* workspace */
    SuiteSparse_long sym    /* 0/1 to indicate whether A is a symmetric matrix */
) ;

/* ============= Routines to get the recommended blocksize ================== */
/* Prototype to get the recommended block size for bidiagonal/tridiagonal
 * reduction of band matrices */
int piro_band_get_blocksize
(
    int m,          /* #rows in the input matrix */
    int n,          /* #columns in the input matrix */
    int kl,         /* lower bandwidth */
    int ku,         /* upper bandwidth */
    int wantuv,     /* flag : 0/1 whether U or V is needed */
    int *blk        /* recommended blocksize, an array of size 4 */
    /* blk[0], blk[1] for #columns, #rows in the block for reducing band in the
     *                upper triangular part.
     * blk[2], blk[3] for #rows, #columns in the block for reducing band in the
     *                lower triangular part.
     **/
) ;

/* Prototype to get the recommended block size for bidiagonal/tridiagonal
 * reduction of band matrices */
int piro_band_get_blocksize_l
(
    SuiteSparse_long m,         /* #rows in the input matrix */
    SuiteSparse_long n,         /* #columns in the input matrix */
    SuiteSparse_long kl,        /* lower bandwidth */
    SuiteSparse_long ku,        /* upper bandwidth */
    SuiteSparse_long wantuv,    /* flag : 0/1 whether U or V is needed */
    SuiteSparse_long *blk       /* recommended blocksize, an array of size 4 */
    /* blk[0], blk[1] for #columns, #rows in the block for reducing band in the
     *                upper triangular part.
     * blk[2], blk[3] for #rows, #columns in the block for reducing band in the
     *                lower triangular part.
     **/
) ;

/* ========================================================================== */
/* === piro_band_check ====================================================== */
/* ========================================================================== */

int piro_band_check_dri     /* return 0 if OK, < 0 if error */
(
    int m,          /* number of rows in A */
    int n,          /* number of columns in A */
    int nrc,        /* number of rows in C */
    int bl,         /* lower bandwidth */
    int bu,         /* upper bandwidth */
    double *A,      /* the band matrix to be reduced */
    int ldab,       /* leading dimension of A */
    double *B1,     /* the diagonal of the bidiagonal matrix */
    double *B2,     /* the superdiagonal of the bidiagional matrix */
    double *U,      /* accumulated left rotations */
    int ldu,        /* leading dimension of U */
    double *V,      /* accumulated right rotations */
    int ldv,        /* leading dimension of V */
    double *C,      /* for apply left rotations */
    int ldc,        /* leading dimension of C */
    int sym,        /* nonzero if A is symmetric, zero if A is unsymmetric */
    int ccode       /* 0: LAPACK interface, nonzero: C interface */
) ;

int piro_band_check_dci
(
    int m,
    int n,
    int nrc,
    int bl,
    int bu,
    double *A,
    int ldab,
    double *B1,
    double *B2,
    double *U,
    int ldu,
    double *V,
    int ldv,
    double *C,
    int ldc,
    int sym,
    int ccode
) ;

int piro_band_check_sri
(
    int m,
    int n,
    int nrc,
    int bl,
    int bu,
    float *A,
    int ldab,
    float *B1,
    float *B2,
    float *U,
    int ldu,
    float *V,
    int ldv,
    float *C,
    int ldc,
    int sym,
    int ccode
) ;

int piro_band_check_sci
(
    int m,
    int n,
    int nrc,
    int bl,
    int bu,
    float *A,
    int ldab,
    float *B1,
    float *B2,
    float *U,
    int ldu,
    float *V,
    int ldv,
    float *C,
    int ldc,
    int sym,
    int ccode
) ;

int piro_band_check_drl
(
    SuiteSparse_long m,
    SuiteSparse_long n,
    SuiteSparse_long nrc,
    SuiteSparse_long bl,
    SuiteSparse_long bu,
    double *A,
    SuiteSparse_long ldab,
    double *B1,
    double *B2,
    double *U,
    SuiteSparse_long ldu,
    double *V,
    SuiteSparse_long ldv,
    double *C,
    SuiteSparse_long ldc,
    SuiteSparse_long sym,
    SuiteSparse_long ccode
) ;

int piro_band_check_dcl
(
    SuiteSparse_long m,
    SuiteSparse_long n,
    SuiteSparse_long nrc,
    SuiteSparse_long bl,
    SuiteSparse_long bu,
    double *A,
    SuiteSparse_long ldab,
    double *B1,
    double *B2,
    double *U,
    SuiteSparse_long ldu,
    double *V,
    SuiteSparse_long ldv,
    double *C,
    SuiteSparse_long ldc,
    SuiteSparse_long sym,
    SuiteSparse_long ccode
) ;

int piro_band_check_srl
(
    SuiteSparse_long m,
    SuiteSparse_long n,
    SuiteSparse_long nrc,
    SuiteSparse_long bl,
    SuiteSparse_long bu,
    float *A,
    SuiteSparse_long ldab,
    float *B1,
    float *B2,
    float *U,
    SuiteSparse_long ldu,
    float *V,
    SuiteSparse_long ldv,
    float *C,
    SuiteSparse_long ldc,
    SuiteSparse_long sym,
    SuiteSparse_long ccode
) ;

int piro_band_check_scl
(
    SuiteSparse_long m,
    SuiteSparse_long n,
    SuiteSparse_long nrc,
    SuiteSparse_long bl,
    SuiteSparse_long bu,
    float *A,
    SuiteSparse_long ldab,
    float *B1,
    float *B2,
    float *U,
    SuiteSparse_long ldu,
    float *V,
    SuiteSparse_long ldv,
    float *C,
    SuiteSparse_long ldc,
    SuiteSparse_long sym,
    SuiteSparse_long ccode
) ;

/* ========================================================================== */
/* ======= Error codes ====================================================== */
/* ========================================================================== */

/* the following error codes match the LAPACK xgbbrd and xxbtrd functions. */
#define PIRO_BAND_OK 0
#define PIRO_BAND_VECT_INVALID (-1)
#define PIRO_BAND_M_INVALID (-2)
#define PIRO_BAND_N_INVALID (-3)
#define PIRO_BAND_NRC_INVALID (-4)
#define PIRO_BAND_BL_INVALID (-5)
#define PIRO_BAND_BU_INVALID (-6)
#define PIRO_BAND_LDAB_INVALID (-8)
#define PIRO_BAND_LDU_INVALID (-12)
#define PIRO_BAND_LDV_INVALID (-14)
#define PIRO_BAND_LDC_INVALID (-16)

/* the following error codes are unique to PIRO_BAND */
#define PIRO_BAND_A_INVALID (-7)
#define PIRO_BAND_B1_INVALID (-9)
#define PIRO_BAND_B2_INVALID (-10)
#define PIRO_BAND_U_INVALID (-11)
#define PIRO_BAND_V_INVALID (-13)
#define PIRO_BAND_C_INVALID (-15)
#define PIRO_BAND_SYM_INVALID (-17)
#define PIRO_BAND_BLKSIZE_INVALID (-18)

#define PIRO_BAND_OUT_OF_MEMORY (-19)
#define PIRO_BAND_UPLO_INVALID (-20)

/* error codes originating from calls to LAPACK */
#define PIRO_BAND_LAPACK_INVALID (-21)
#define PIRO_BAND_LAPACK_FAILURE (-22)

/* ========================================================================== */
/* === flop count =========================================================== */
/* ========================================================================== */

/* The flop count is not computed unless PIRO_BAND is compiled with
   -DBENCHMARK.  The flop count assumes the matrix is real.  If the matrix
   is complex, the flop count is roughly 4 times higher. */

/* TODO rename BENCHMARK to PIRO_BAND_BENCHMARK */
#ifdef BENCHMARK
extern double piro_band_flops ;
#define PIRO_BAND_FLOPS(fl)   { piro_band_flops += ( ((fl) > 0) ? (fl): 0) ; }
#define PIRO_BAND_CLEAR_FLOPS { piro_band_flops = 0 ; }
#define PIRO_BAND_GET_FLOPS   (piro_band_flops)
#else
#define PIRO_BAND_FLOPS(fl)   ;
#define PIRO_BAND_CLEAR_FLOPS ;
#define PIRO_BAND_GET_FLOPS   (0)
#endif /* BENCHMARK */

/* ========================================================================== */
/* === version ============================================================== */
/* ========================================================================== */

/* TODO pick a date */
#define PIRO_BAND_DATE "May 9, 2014"
#define PIRO_BAND_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define PIRO_BAND_MAIN_VERSION 1
#define PIRO_BAND_SUB_VERSION 0
#define PIRO_BAND_SUBSUB_VERSION 0
#define PIRO_BAND_VERSION PIRO_BAND_VERSION_CODE(PIRO_BAND_MAIN_VERSION,PIRO_BAND_SUB_VERSION)

#ifdef __cplusplus
}
#endif

#endif /* PIRO_BAND_H */
