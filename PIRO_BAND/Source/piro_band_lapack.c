/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_lapack.c ================================== */
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
 * Replacement functions for LAPACK's band reduction. Prefix all LAPACK names
 * with "piro_band_" to get the piro_band equivalent of the LAPACK function.
 *
 * For example, the piro_band_dgbbrd function is a replacement for the LAPACK
 * function dgbbrd. The arguments and the return values are the same as in
 * LAPACK. The long versions of LAPACK functions will have a suffix _l like
 * piro_band_dgbbrd_l.  The LAPACK functions whose equivalents are in this
 * file are:
 *
 * DGBBRD: reduce a real band matrix to upper bidiagonal form
 * SGBBRD: reduce a real band matrix to upper bidiagonal form
 * ZGBBRD: reduce a complex band matrix to real upper bidiagonal form
 * CGBBRD: reduce a complex band matrix to real upper bidiagonal form
 * DSBTRD: reduce a real sym. band matrix to sym. tridiagonal
 * SSBTRD: reduce a real sym. band matrix to sym. tridiagonal
 * ZHBTRD: reduce a complex Hermitian band matrix to real sym. tridiagonal
 * CHBTRD: reduce a complex Hermitian band matrix to real sym. tridiagonal
 * */

#include "piro_band_memory.h"
#include "piro_band_lapack.h"
#include "piro_band_lapack_internal.h"
#include "piro_band.h"

/* ========================================================================== */
/* DGBBRD: reduce a real band matrix to upper bidiagonal form */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(dgbbrd)
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  Int m,            /* #rows in AB */
  Int n,            /* #columns in AB */
  Int ncc,          /* #columns in C */
  Int kl,           /* lower bandwidth */
  Int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output: super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated left rotations */
  Int ldq,          /* leading dimension of Q */
  double *PT,       /* Accumulated right rotations */
  Int ldpt,         /* leading dimension of PT */
  double *C,        /* Accumulated left rotations */
  Int ldc,          /* leading dimension of C */
  double *work,     /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int wantq, wantpt, wantb, msize ;
    Int blk[4] ;
    double *dwork ;
    double *A, *newA, *newC ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    wantb = (vect[0] == 'B') ;                  /* both Q and PT */
    wantq = ((vect[0] == 'Q') || wantb) ;       /* just Q */
    wantpt = ((vect[0] == 'P') || wantb) ;      /* just PT */
    if (!wantq && !wantpt && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (wantpt && PT == NULL) *info = PIRO_BAND_V_INVALID ;
    if (ncc > 0 && C == NULL) *info = PIRO_BAND_C_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_dr) (m, n, ncc, kl, ku, AB, ldab,
        D, E, Q, ldq, PT, ldpt, C, ldc, 0, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    /* transpose C if required */
    newC = NULL ;
    if (C != NULL && ncc > 0)
    {
        newC = (double *) MALLOC (MAX (1, ncc * m) * sizeof(double)) ;
        if (newC == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(general_transpose_dr) (m, ncc, C, ldc, newC, ncc);
    }

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (double *) MALLOC(MAX (1, n * (kl + 2)) * sizeof(double)) ;
        if (newA == NULL)
        {
            if (newC != NULL)
            {
                FREE(newC) ;
            }
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_dr) (n, AB, ldab, kl, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = kl + 2 ;
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (m, n, kl, ku,
                        (wantq || wantpt || ncc > 0), blk) ;

    /* Allocate work for the band reduction*/
    msize = 2 * MAX (blk[0]*blk[1], blk[2]*blk[3]) ;
    dwork = (double *) MALLOC (MAX (1,msize) * sizeof(double)) ;

    if (dwork == NULL && (kl != 0 || ku > 1))
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        /* Reduce the matrix to bidiagonal */
        *info = PIRO_BAND_LAPACK_NAME(reduce_dr) (blk, m, n, ncc, kl, ku, A,
            ldab, D, E, Q, ldq, PT, ldpt, newC, ncc, dwork, 0) ;

        if (wantpt)
        {
            /* Transpose PT as reduce finds P */
            PIRO_BAND_LAPACK_NAME(inplace_conjugate_transpose_dr) (n, PT, ldpt);
        }

        if (C != NULL && ncc > 0)
        {
            /* Transpose C back */
            PIRO_BAND_LAPACK_NAME(general_transpose_dr) (ncc, m, newC, ncc, C,
                    ldc) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newC != NULL)
    {
        FREE(newC) ;
    }
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}

/* ========================================================================== */
/* SGBBRD: reduce a real band matrix to upper bidiagonal form */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(sgbbrd)
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  Int m,            /* #rows in AB */
  Int n,            /* #columns in AB */
  Int ncc,          /* #columns in C */
  Int kl,           /* lower bandwidth */
  Int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output: super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated left rotations */
  Int ldq,          /* leading dimension of Q */
  float *PT,        /* Accumulated right rotations */
  Int ldpt,         /* leading dimension of PT */
  float *C,         /* Accumulated left rotations */
  Int ldc,          /* leading dimension of C */
  float *work,      /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int wantq, wantpt, wantb, msize ;
    Int blk[4] ;
    float *dwork ;
    float *A, *newA, *newC ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    wantb = (vect[0] == 'B') ;                  /* both Q and PT */
    wantq = ((vect[0] == 'Q') || wantb) ;       /* just Q */
    wantpt = ((vect[0] == 'P') || wantb) ;      /* just PT */
    if (!wantq && !wantpt && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (wantpt && PT == NULL) *info = PIRO_BAND_V_INVALID ;
    if (ncc > 0 && C == NULL) *info = PIRO_BAND_C_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_sr) (m, n, ncc, kl, ku, AB, ldab,
        D, E, Q, ldq, PT, ldpt, C, ldc, 0, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    /* transpose C if required */
    newC = NULL ;
    if (C != NULL && ncc > 0)
    {
        newC = (float *) MALLOC(MAX (1, ncc * m) * sizeof(float)) ;
        if (newC == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(general_transpose_sr) (m, ncc, C, ldc, newC, ncc);
    }

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (float *) MALLOC(MAX (1, n * (kl + 2)) * sizeof(float)) ;
        if (newA == NULL)
        {
            if (newC != NULL)
            {
                FREE(newC) ;
            }
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_sr) (n, AB, ldab, kl, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = kl + 2 ;
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (m, n, kl, ku,
                        (wantq || wantpt || ncc > 0), blk) ;

    /* Allocate work for the band reduction*/
    msize = 2 * MAX (blk[0]*blk[1], blk[2]*blk[3]) ;
    dwork = (float *) MALLOC (MAX (1,msize) * sizeof(float)) ;

    if (dwork == NULL && (kl != 0 || ku > 1))
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        /* Reduce the matrix to bidiagonal */
        *info = PIRO_BAND_LAPACK_NAME(reduce_sr) (blk, m, n, ncc, kl, ku, A,
            ldab, D, E, Q, ldq, PT, ldpt, newC, ncc, dwork, 0) ;

        if (wantpt)
        {
            /* Transpose PT as reduce finds P */
            PIRO_BAND_LAPACK_NAME(inplace_conjugate_transpose_sr) (n, PT, ldpt);
        }

        if (C != NULL && ncc > 0)
        {
            /* Transpose C back */
            PIRO_BAND_LAPACK_NAME(general_transpose_sr) (ncc, m, newC, ncc, C,
                    ldc) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newC != NULL)
    {
        FREE(newC) ;
    }
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}

/* ========================================================================== */
/* ZGBBRD: reduce a complex band matrix to real upper bidiagonal form */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(zgbbrd)
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  Int m,            /* #rows in AB */
  Int n,            /* #columns in AB */
  Int ncc,          /* #columns in C */
  Int kl,           /* lower bandwidth */
  Int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output: super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated left rotations */
  Int ldq,          /* leading dimension of Q */
  double *PT,       /* Accumulated right rotations */
  Int ldpt,         /* leading dimension of PT */
  double *C,        /* Accumulated left rotations */
  Int ldc,          /* leading dimension of C */
  double *work,     /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int wantq, wantpt, wantb, msize ;
    Int blk[4] ;
    double *dwork ;
    double *A, *newA, *newC ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    wantb = (vect[0] == 'B') ;                  /* both Q and PT */
    wantq = ((vect[0] == 'Q') || wantb) ;       /* just Q */
    wantpt = ((vect[0] == 'P') || wantb) ;      /* just PT */
    if (!wantq && !wantpt && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (wantpt && PT == NULL) *info = PIRO_BAND_V_INVALID ;
    if (ncc > 0 && C == NULL) *info = PIRO_BAND_C_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_dc) (m, n, ncc, kl, ku, AB, ldab,
        D, E, Q, ldq, PT, ldpt, C, ldc, 0, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    /* transpose C if required */
    newC = NULL ;
    if (C != NULL && ncc > 0)
    {
        newC = (double *) MALLOC(MAX (1, 2 * ncc * m) * sizeof(double)) ;
        if (newC == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(general_transpose_dc) (m, ncc, C, ldc, newC, ncc);
    }

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (double *) MALLOC (MAX (1, 2 * n * (kl + 2)) * sizeof(double)) ;
        if (newA == NULL)
        {
            if (newC != NULL)
            {
                FREE(newC) ;
            }
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_dc) (n, AB, ldab, kl, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = kl + 2 ;
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (m, n, kl, ku,
                        (wantq || wantpt || ncc > 0), blk) ;

    /* Allocate work for the band reduction*/
    msize = 4 * MAX (blk[0]*blk[1], blk[2]*blk[3]) ;
    dwork = (double *) MALLOC (MAX (1,msize) * sizeof(double)) ;

    if (dwork == NULL && (kl != 0 || ku > 1))
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        /* Reduce the matrix to bidiagonal */
        *info = PIRO_BAND_LAPACK_NAME(reduce_dc) (blk, m, n, ncc, kl, ku, A,
            ldab, D, E, Q, ldq, PT, ldpt, newC, ncc, dwork, 0) ;

        if (wantpt)
        {
            /* Transpose PT as reduce finds P */
            PIRO_BAND_LAPACK_NAME(inplace_conjugate_transpose_dc) (n, PT, ldpt);
        }

        if (C != NULL && ncc > 0)
        {
            /* Transpose C back */
            PIRO_BAND_LAPACK_NAME(general_transpose_dc) (ncc, m, newC, ncc, C,
                    ldc) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newC != NULL)
    {
        FREE(newC) ;
    }
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}

/* ========================================================================== */
/* CGBBRD: reduce a complex band matrix to real upper bidiagonal form */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(cgbbrd)
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  Int m,            /* #rows in AB */
  Int n,            /* #columns in AB */
  Int ncc,          /* #columns in C */
  Int kl,           /* lower bandwidth */
  Int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output: super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated left rotations */
  Int ldq,          /* leading dimension of Q */
  float *PT,        /* Accumulated right rotations */
  Int ldpt,         /* leading dimension of PT */
  float *C,         /* Accumulated left rotations */
  Int ldc,          /* leading dimension of C */
  float *work,      /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int wantq, wantpt, wantb, msize ;
    Int blk[4] ;
    float *dwork ;
    float *A, *newA, *newC ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    wantb = (vect[0] == 'B') ;                  /* both Q and PT */
    wantq = ((vect[0] == 'Q') || wantb) ;       /* just Q */
    wantpt = ((vect[0] == 'P') || wantb) ;      /* just PT */
    if (!wantq && !wantpt && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (wantpt && PT == NULL) *info = PIRO_BAND_V_INVALID ;
    if (ncc > 0 && C == NULL) *info = PIRO_BAND_C_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_sc) (m, n, ncc, kl, ku, AB, ldab,
        D, E, Q, ldq, PT, ldpt, C, ldc, 0, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    /* transpose C if required */
    newC = NULL ;
    if (C != NULL && ncc > 0)
    {
        newC = (float *) MALLOC(MAX (1, 2 * ncc * m) * sizeof(float)) ;
        if (newC == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(general_transpose_sc) (m, ncc, C, ldc, newC, ncc);
    }

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (float *) MALLOC(MAX (1, 2 * n * (kl + 2)) * sizeof(float)) ;
        if (newA == NULL)
        {
            if (newC != NULL)
            {
                FREE(newC) ;
            }
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_sc) (n, AB, ldab, kl, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = kl + 2 ;
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (m, n, kl, ku,
                        (wantq || wantpt || ncc > 0), blk) ;

    /* Allocate work for the band reduction*/
    msize = 4 * MAX (blk[0]*blk[1], blk[2]*blk[3]) ;
    dwork = (float *) MALLOC (MAX (1,msize) * sizeof(float)) ;

    if (dwork == NULL && (kl != 0 || ku > 1))
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        /* Reduce the matrix to bidiagonal */
        *info = PIRO_BAND_LAPACK_NAME(reduce_sc) (blk, m, n, ncc, kl, ku, A,
            ldab, D, E, Q, ldq, PT, ldpt, newC, ncc, dwork, 0) ;

        if (wantpt)
        {
            /* Transpose PT as reduce finds P */
            PIRO_BAND_LAPACK_NAME(inplace_conjugate_transpose_sc) (n, PT, ldpt);
        }

        if (C != NULL && ncc > 0)
        {
            /* Transpose C back */
            PIRO_BAND_LAPACK_NAME(general_transpose_sc) (ncc, m, newC, ncc, C,
                    ldc) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newC != NULL)
    {
        FREE(newC) ;
    }
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}

/* ========================================================================== */
/* DSBTRD: reduce a real sym. band matrix to sym. tridiagonal */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(dsbtrd)
(
  char *vect,       /* 'N': do not form Q, 'V': form Q, 'U': update Q */
  char *uplo,       /* 'U': upper triangle of A stored, 'L': lower */
  Int n,            /* #columns in AB */
  Int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output : super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated rotations */
  Int ldq,          /* leading dimension of Q */
  double *work,     /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int wantq, initq ;
    Int blk[4] ;
    double *dwork, *upper ;
    Int ncc ;
    double *A, *newA ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    initq = (vect[0] == 'V') ;
    wantq = ((vect[0] == 'U') || initq) ;
    if (!wantq && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (uplo[0] != 'U' && uplo[0] != 'L') *info = PIRO_BAND_UPLO_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_dr) (n, n, 0, 0, ku, AB, ldab,
        D, E, Q, ldq, NULL, 0, NULL, 0, 1, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    ncc = wantq ? n : 0 ;

    /* Initialize Q to identity */
    if (initq) PIRO_BAND_LAPACK_NAME(set_to_eye_dr) (ldq, n, Q) ;

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (double *) MALLOC(MAX (1, n * 2) * sizeof(double)) ;
        if (newA == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_dr) (n, AB, ldab, ku, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = 2 ;
        uplo[0] = 'U' ; /* We stored the upper diagonal */
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (n, n, 0, ku, wantq, blk) ;

    /* Allocate work for the band reduction*/
    dwork = (double *) MALLOC(MAX (1, 2 * blk[0] * blk[1]) * sizeof(double)) ;
    if (dwork == NULL && ku > 1)
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        if (uplo[0] == 'L')
        {
            /* If lower triangular part of the band is stores transpose it */
            upper = (double *) MALLOC(MAX (1, n * (ku+1)) * sizeof(double)) ;
            if (upper == NULL)
            {
                *info = PIRO_BAND_OUT_OF_MEMORY ;
                FREE (dwork) ;
                ASSERT (newA == NULL) ;
                return ;
            }
            PIRO_BAND_LAPACK_NAME(lowerband_transpose_dr) (n, ku, A, ldab,
                upper) ;

            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_dr) (blk, n, n, ncc, 0, ku,
                upper, ku+1, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;

            FREE(upper) ;
        }
        else
        {
            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_dr) (blk, n, n, ncc, 0, ku, A,
                    ldab, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}

/* ========================================================================== */
/* SSBTRD: reduce a real sym. band matrix to sym. tridiagonal */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(ssbtrd)
(
  char *vect,       /* 'N': do not form Q, 'V': form Q, 'U': update Q */
  char *uplo,       /* 'U': upper triangle of A stored, 'L': lower */
  Int n,            /* #columns in AB */
  Int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output : super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated rotations */
  Int ldq,          /* leading dimension of Q */
  float *work,      /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int initq, wantq ;
    Int blk[4] ;
    float *dwork, *upper ;
    Int ncc ;
    float *A, *newA ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    initq = (vect[0] == 'V') ;
    wantq = ((vect[0] == 'U') || initq) ;
    if (!wantq && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (uplo[0] != 'U' && uplo[0] != 'L') *info = PIRO_BAND_UPLO_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_sr) (n, n, 0, 0, ku, AB, ldab,
        D, E, Q, ldq, NULL, 0, NULL, 0, 1, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    ncc = wantq ? n : 0 ;

    /* Initialize Q to identity */
    if (initq) PIRO_BAND_LAPACK_NAME(set_to_eye_sr) (ldq, n, Q) ;

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (float *) MALLOC(MAX (1, n * 2) * sizeof(float)) ;
        if (newA == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_sr) (n, AB, ldab, ku, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = 2 ;
        uplo[0] = 'U' ; /* We stored the upper diagonal */
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (n, n, 0, ku, wantq, blk) ;

    /* Allocate work for the band reduction*/
    dwork = (float *) MALLOC(MAX (1, 2 * blk[0] * blk[1]) * sizeof(float)) ;
    if (dwork == NULL && ku > 1)
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        if (uplo[0] == 'L')
        {
            /* If lower triangular part of the band is stores transpose it */
            upper = (float *) MALLOC(MAX (1, n * (ku+1)) * sizeof(float)) ;
            if (upper == NULL)
            {
                *info = PIRO_BAND_OUT_OF_MEMORY ;
                FREE (dwork) ;
                ASSERT (newA == NULL) ;
                return ;
            }
            PIRO_BAND_LAPACK_NAME(lowerband_transpose_sr) (n, ku, A, ldab,
                upper) ;

            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_sr) (blk, n, n, ncc, 0, ku,
                upper, ku+1, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;

            FREE(upper) ;
        }
        else
        {
            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_sr) (blk, n, n, ncc, 0, ku, A,
                    ldab, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}

/* ========================================================================== */
/* ZHBTRD: reduce a complex Hermitian band matrix to real sym. tridiagonal */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(zhbtrd)
(
  char *vect,       /* 'N': do not form Q, 'V': form Q, 'U': update Q */
  char *uplo,       /* 'U': upper triangle of A stored, 'L': lower */
  Int n,            /* #columns in AB */
  Int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output : super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated rotations */
  Int ldq,          /* leading dimension of Q */
  double *work,     /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int initq, wantq ;
    Int blk[4] ;
    double *dwork, *upper ;
    Int ncc ;
    double *A, *newA ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    initq = (vect[0] == 'V') ;
    wantq = ((vect[0] == 'U') || initq) ;
    if (!wantq && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (uplo[0] != 'U' && uplo[0] != 'L') *info = PIRO_BAND_UPLO_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_dc) (n, n, 0, 0, ku, AB, ldab,
        D, E, Q, ldq, NULL, 0, NULL, 0, 1, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    ncc = wantq ? n : 0 ;

    /* Initialize Q to identity */
    if (initq) PIRO_BAND_LAPACK_NAME(set_to_eye_dc) (ldq, n, Q) ;

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (double *) MALLOC(MAX (1, 2 * n * 2) * sizeof(double)) ;
        if (newA == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_dc) (n, AB, ldab, ku, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = 2 ;
        uplo[0] = 'U' ; /* We stored the upper diagonal */
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (n, n, 0, ku, wantq, blk) ;

    /* Allocate work for the band reduction*/
    dwork = (double *) MALLOC(MAX (1, 2*2 * blk[0] * blk[1]) * sizeof(double)) ;
    if (dwork == NULL && ku > 1)
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        if (uplo[0] == 'L')
        {
            /* If lower triangular part of the band is stores transpose it */
            upper = (double *)MALLOC(MAX (1, 2 * n * (ku+1)) * sizeof(double));
            if (upper == NULL)
            {
                *info = PIRO_BAND_OUT_OF_MEMORY ;
                FREE (dwork) ;
                ASSERT (newA == NULL) ;
                return ;
            }
            PIRO_BAND_LAPACK_NAME(lowerband_transpose_dc) (n, ku, A, ldab,
                upper) ;

            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_dc) (blk, n, n, ncc, 0, ku,
                upper, ku+1, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;

            FREE(upper) ;
        }
        else
        {
            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_dc) (blk, n, n, ncc, 0, ku, A,
                    ldab, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}

/* ========================================================================== */
/* CHBTRD: reduce a complex Hermitian band matrix to real sym. tridiagonal */
/* ========================================================================== */

void PIRO_BAND_LONG_NAME(chbtrd)
(
  char *vect,       /* 'N': do not form Q, 'V': form Q, 'U': update Q */
  char *uplo,       /* 'U': upper triangle of A stored, 'L': lower */
  Int n,            /* #columns in AB */
  Int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  Int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output : super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated rotations */
  Int ldq,          /* leading dimension of Q */
  float *work,      /* workspace : ignored now */
  Int *info         /* info = 0 for success and negative for failure */
)
{
    Int initq, wantq ;
    Int blk[4] ;
    float *dwork, *upper ;
    Int ncc ;
    float *A, *newA ;

    /* check inputs */
    *info = PIRO_BAND_OK ;
    initq = (vect[0] == 'V') ;
    wantq = ((vect[0] == 'U') || initq) ;
    if (!wantq && vect[0] != 'N') *info = PIRO_BAND_VECT_INVALID ;
    if (uplo[0] != 'U' && uplo[0] != 'L') *info = PIRO_BAND_UPLO_INVALID ;
    if (wantq  && Q  == NULL) *info = PIRO_BAND_U_INVALID ;
    if (*info != PIRO_BAND_OK) return ;
    *info = PIRO_BAND_LAPACK_NAME (check_sc) (n, n, 0, 0, ku, AB, ldab,
        D, E, Q, ldq, NULL, 0, NULL, 0, 1, 0) ;
    if (*info != PIRO_BAND_OK) return ;

    ncc = wantq ? n : 0 ;

    /* Initialize Q to identity */
    if (initq) PIRO_BAND_LAPACK_NAME(set_to_eye_sc) (ldq, n, Q) ;

    A = AB ;
    newA = NULL ;
    if (ku == 0)
    {
        /* Add a upper bidiagonal with zeros to AB */
        newA = (float *) MALLOC(MAX (1, 2 * n * 2) * sizeof(float)) ;
        if (newA == NULL)
        {
            *info = PIRO_BAND_OUT_OF_MEMORY ;
            return ;
        }
        PIRO_BAND_LAPACK_NAME(add_bidiag_sc) (n, AB, ldab, ku, newA) ;
        A = newA ;
        ku = 1 ;
        ldab = 2 ;
        uplo[0] = 'U' ; /* We stored the upper diagonal */
    }

    /* Find block size */
    PIRO_BAND_LONG_NAME(get_blocksize) (n, n, 0, ku, wantq, blk) ;

    /* Allocate work for the band reduction*/
    dwork = (float *) MALLOC(MAX (1, 2*2 * blk[0] * blk[1]) * sizeof(float)) ;
    if (dwork == NULL && ku > 1)
    {
        *info = PIRO_BAND_OUT_OF_MEMORY ;
    }
    else
    {
        if (uplo[0] == 'L')
        {
            /* If lower triangular part of the band is stores transpose it */
            upper = (float *) MALLOC(MAX (1, 2 * n * (ku+1)) * sizeof(float)) ;
            if (upper == NULL)
            {
                *info = PIRO_BAND_OUT_OF_MEMORY ;
                FREE (dwork) ;
                ASSERT (newA == NULL) ;
                return ;
            }
            PIRO_BAND_LAPACK_NAME(lowerband_transpose_sc) (n, ku, A, ldab,
                upper) ;

            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_sc) (blk, n, n, ncc, 0, ku,
                upper, ku+1, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;

            FREE(upper) ;
        }
        else
        {
            /* Reduce the matrix to tridiagonal */
            *info = PIRO_BAND_LAPACK_NAME(reduce_sc) (blk, n, n, ncc, 0, ku, A,
                    ldab, D, E, NULL, 0, NULL, 0, Q, ldq, dwork, 1) ;
        }
        FREE(dwork) ;
    }

    /* free work */
    if (newA != NULL)
    {
        FREE(newA) ;
    }

    return ;
}
