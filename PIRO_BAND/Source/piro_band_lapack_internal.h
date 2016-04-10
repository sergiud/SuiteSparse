/* ========================================================================== */
/* === PIRO_BAND/Include/piro_band_lapack_internal.h ======================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* This file should not be #include'd in user programs. */

#ifndef PIRO_BAND_LAPACK_INTERNAL_H
#define PIRO_BAND_LAPACK_INTERNAL_H


/* ========================================================================== */
/* ===================== Definitions for internal routines  ================= */
/* ========================================================================== */

void piro_band_set_to_eye_sri
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    float *A        /* Matrix to be initialized to identity */
) ;

void piro_band_set_to_eye_sci
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    float *A        /* Matrix to be initialized to identity */
) ;

void piro_band_set_to_eye_dri
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    double *A       /* Matrix to be initialized to identity */
) ;

void piro_band_set_to_eye_dci
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    double *A       /* Matrix to be initialized to identity */
) ;

void piro_band_set_to_eye_srl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    float *A                /* Matrix to be initialized to identity */
) ;

void piro_band_set_to_eye_scl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    float *A                /* Matrix to be initialized to identity */
) ;

void piro_band_set_to_eye_drl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    double *A               /* Matrix to be initialized to identity */
) ;

void piro_band_set_to_eye_dcl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    double *A               /* Matrix to be initialized to identity */
) ;

void piro_band_general_transpose_dri
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    double *A,      /* Matrix to be transposed */
    int lda,        /* leading dimension of A */
    double *AT,     /* A' */
    int ldat        /* leading dimension of A */
) ;

void piro_band_general_transpose_dci
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    double *A,      /* Matrix to be transposed */
    int lda,        /* leading dimension of A */
    double *AT,     /* A' */
    int ldat        /* leading dimension of A */
) ;

void piro_band_general_transpose_sri
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    float *A,       /* Matrix to be transposed */
    int lda,        /* leading dimension of A */
    float *AT,      /* A' */
    int ldat        /* leading dimension of A */
) ;

void piro_band_general_transpose_sci
(
    int m,          /* #rows in A */
    int n,          /* #columns in A */
    float *A,       /* Matrix to be transposed */
    int lda,        /* leading dimension of A */
    float *AT,      /* A' */
    int ldat        /* leading dimension of A */
) ;

void piro_band_general_transpose_drl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    double *A,              /* Matrix to be transposed */
    SuiteSparse_long lda,   /* leading dimension of A */
    double *AT,             /* A' */
    SuiteSparse_long ldat   /* leading dimension of A */
) ;

void piro_band_general_transpose_dcl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    double *A,              /* Matrix to be transposed */
    SuiteSparse_long lda,   /* leading dimension of A */
    double *AT,             /* A' */
    SuiteSparse_long ldat   /* leading dimension of A */
) ;

void piro_band_general_transpose_srl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    float *A,               /* Matrix to be transposed */
    SuiteSparse_long lda,   /* leading dimension of A */
    float *AT,              /* A' */
    SuiteSparse_long ldat   /* leading dimension of A */
) ;

void piro_band_general_transpose_scl
(
    SuiteSparse_long m,     /* #rows in A */
    SuiteSparse_long n,     /* #columns in A */
    float *A,               /* Matrix to be transposed */
    SuiteSparse_long lda,   /* leading dimension of A */
    float *AT,              /* A' */
    SuiteSparse_long ldat   /* leading dimension of A */
) ;

void piro_band_lowerband_transpose_dri
(
    int n,          /* #columns */
    int bl,         /* lower bandwidth */
    double *A,      /* band matrix in packed band format */
    int lda,        /* leading dimension of A */
    double *AT      /* A' in upper band format */
) ;

void piro_band_lowerband_transpose_dci
(
    int n,          /* #columns */
    int bl,         /* lower bandwidth */
    double *A,      /* band matrix in packed band format */
    int lda,        /* leading dimension of A */
    double *AT      /* A' in upper band format */
) ;

void piro_band_lowerband_transpose_sri
(
    int n,          /* #columns */
    int bl,         /* lower bandwidth */
    float *A,       /* band matrix in packed band format */
    int lda,        /* leading dimension of A */
    float *AT       /* A' in upper band format */
) ;

void piro_band_lowerband_transpose_sci
(
    int n,          /* #columns */
    int bl,         /* lower bandwidth */
    float *A,       /* band matrix in packed band format */
    int lda,        /* leading dimension of A */
    float *AT       /* A' in upper band format */
) ;

void piro_band_lowerband_transpose_drl
(
    SuiteSparse_long n,     /* #columns */
    SuiteSparse_long bl,    /* lower bandwidth */
    double *A,              /* band matrix in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    double *AT              /* A' in upper band format */
) ;

void piro_band_lowerband_transpose_dcl
(
    SuiteSparse_long n,     /* #columns */
    SuiteSparse_long bl,    /* lower bandwidth */
    double *A,              /* band matrix in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    double *AT              /* A' in upper band format */
) ;

void piro_band_lowerband_transpose_srl
(
    SuiteSparse_long n,     /* #columns */
    SuiteSparse_long bl,    /* lower bandwidth */
    float *A,               /* band matrix in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    float *AT               /* A' in upper band format */
) ;

void piro_band_lowerband_transpose_scl
(
    SuiteSparse_long n,     /* #columns */
    SuiteSparse_long bl,    /* lower bandwidth */
    float *A,               /* band matrix in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    float *AT               /* A' in upper band format */
) ;

void piro_band_shift_superdiag_dri
(
    int n,          /* #columns */
    double *E       /* first super diagonal */
) ;

void piro_band_shift_superdiag_dci
(
    int n,          /* #columns */
    double *E       /* first super diagonal */
) ;

void piro_band_shift_superdiag_sri
(
    int n,          /* #columns */
    float *E        /* first super diagonal */
) ;

void piro_band_shift_superdiag_sci
(
    int n,          /* #columns */
    float *E        /* first super diagonal */
) ;

void piro_band_shift_superdiag_drl
(
    SuiteSparse_long n,     /* #columns */
    double *E               /* first super diagonal */
) ;

void piro_band_shift_superdiag_dcl
(
    SuiteSparse_long n,     /* #columns */
    double *E               /* first super diagonal */
) ;

void piro_band_shift_superdiag_srl
(
    SuiteSparse_long n,     /* #columns */
    float *E                /* first super diagonal */
) ;

void piro_band_shift_superdiag_scl
(
    SuiteSparse_long n,     /* #columns */
    float *E                /* first super diagonal */
) ;

void piro_band_add_bidiag_dri
(
    int n,          /* #columns */
    double *A,      /* A stored in packed band format */
    int lda,        /* leading dimension of A */
    int kl,         /* lower bandwidth of A */
    double *newA    /* A with the diagonal and first super diagonal */
) ;

void piro_band_add_bidiag_dci
(
    int n,          /* #columns */
    double *A,      /* A stored in packed band format */
    int lda,        /* leading dimension of A */
    int kl,         /* lower bandwidth of A */
    double *newA    /* A with the diagonal and first super diagonal */
) ;

void piro_band_add_bidiag_sri
(
    int n,          /* #columns */
    float *A,       /* A stored in packed band format */
    int lda,        /* leading dimension of A */
    int kl,         /* lower bandwidth of A */
    float *newA     /* A with the diagonal and first super diagonal */
) ;

void piro_band_add_bidiag_sci
(
    int n,          /* #columns */
    float *A,       /* A stored in packed band format */
    int lda,        /* leading dimension of A */
    int kl,         /* lower bandwidth of A */
    float *newA     /* A with the diagonal and first super diagonal */
) ;

void piro_band_add_bidiag_drl
(
    SuiteSparse_long n,     /* #columns */
    double *A,              /* A stored in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    SuiteSparse_long kl,    /* lower bandwidth of A */
    double *newA            /* A with the diagonal and first super diagonal */
) ;

void piro_band_add_bidiag_dcl
(
    SuiteSparse_long n,     /* #columns */
    double *A,              /* A stored in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    SuiteSparse_long kl,    /* lower bandwidth of A */
    double *newA            /* A with the diagonal and first super diagonal */
) ;

void piro_band_add_bidiag_srl
(
    SuiteSparse_long n,     /* #columns */
    float *A,               /* A stored in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    SuiteSparse_long kl,    /* lower bandwidth of A */
    float *newA             /* A with the diagonal and first super diagonal */
) ;

void piro_band_add_bidiag_scl
(
    SuiteSparse_long n,     /* #columns */
    float *A,               /* A stored in packed band format */
    SuiteSparse_long lda,   /* leading dimension of A */
    SuiteSparse_long kl,    /* lower bandwidth of A */
    float *newA             /* A with the diagonal and first super diagonal */
) ;

void piro_band_inplace_conjugate_transpose_dri
(
    int n,          /* #columns */
    double *A,      /* A stored in packed band format */
    int lda         /* leading dimension of A */
) ;

void piro_band_inplace_conjugate_transpose_dci
(
    int n,          /* #columns */
    double *A,      /* A stored in packed band format */
    int lda         /* leading dimension of A */
) ;

void piro_band_inplace_conjugate_transpose_sri
(
    int n,          /* #columns */
    float *A,       /* A stored in packed band format */
    int lda         /* leading dimension of A */
) ;

void piro_band_inplace_conjugate_transpose_sci
(
    int n,          /* #columns */
    float *A,       /* A stored in packed band format */
    int lda         /* leading dimension of A */
) ;

void piro_band_inplace_conjugate_transpose_drl
(
    SuiteSparse_long n,     /* #columns */
    double *A,              /* A stored in packed band format */
    SuiteSparse_long lda    /* leading dimension of A */
) ;

void piro_band_inplace_conjugate_transpose_dcl
(
    SuiteSparse_long n,     /* #columns */
    double *A,              /* A stored in packed band format */
    SuiteSparse_long lda    /* leading dimension of A */
) ;

void piro_band_inplace_conjugate_transpose_srl
(
    SuiteSparse_long n,     /* #columns */
    float *A,               /* A stored in packed band format */
    SuiteSparse_long lda    /* leading dimension of A */
) ;

void piro_band_inplace_conjugate_transpose_scl
(
    SuiteSparse_long n,     /* #columns */
    float *A,               /* A stored in packed band format */
    SuiteSparse_long lda    /* leading dimension of A */
) ;

#endif /* PIRO_BAND_LAPACK_INTERNAL_H */
