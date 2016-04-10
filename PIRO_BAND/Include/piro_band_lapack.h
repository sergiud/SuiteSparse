/* ========================================================================== */
/* === PIRO_BAND/Include/piro_band_lapack.h ================================= */
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
 * External interface for replacement functions for LAPACK's band reduction
 * routines.  To obtain the function name prefix all LAPACK routine names
 * with "piro_band_" to get the piro_band equivalent of it.
 *
 * For example, dgbbrd will piro_band_dgbbrd(). The arguments and the return
 * values are the same as in LAPACK. The long versions of LAPACK functions will
 * have a suffix _l like piro_band_dgbbrd_l().
 *
 * */

#ifndef PIRO_BAND_LAPACK_H
#define PIRO_BAND_LAPACK_H

/* make it easy for C++ programs to include piro_band_reduce */
#ifdef __cplusplus
extern "C" {
#endif

#include "SuiteSparse_config.h"

void piro_band_dgbbrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  int m,            /* #rows in AB */
  int n,            /* #columns in AB */
  int ncc,          /* #columns in C */
  int kl,           /* lower bandwidth */
  int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output : super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated left rotations */
  int ldq,          /* leading dimension of Q */
  double *PT,       /* Accumulated right rotations */
  int ldpt,         /* leading dimension of PT */
  double *C,        /* Accumulated left rotations */
  int ldc,          /* leading dimension of C */
  double *work,     /* workspace : ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;

void piro_band_sgbbrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  int m,            /* #rows in AB */
  int n,            /* #columns in AB */
  int ncc,          /* #columns in C */
  int kl,           /* lower bandwidth */
  int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output : super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated left rotations */
  int ldq,          /* leading dimension of Q */
  float *PT,        /* Accumulated right rotations */
  int ldpt,         /* leading dimension of PT */
  float *C,         /* Accumulated left rotations */
  int ldc,          /* leading dimension of C */
  float *work,      /* workspace : ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;

void piro_band_zgbbrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  int m,            /* #rows in AB */
  int n,            /* #columns in AB */
  int ncc,          /* #columns in C */
  int kl,           /* lower bandwidth */
  int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output : super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated left rotations */
  int ldq,          /* leading dimension of Q */
  double *PT,       /* Accumulated right rotations */
  int ldpt,         /* leading dimension of PT */
  double *C,        /* Accumulated left rotations */
  int ldc,          /* leading dimension of C */
  double *work,     /* workspace : ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;

void piro_band_cgbbrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  int m,            /* #rows in AB */
  int n,            /* #columns in AB */
  int ncc,          /* #columns in C */
  int kl,           /* lower bandwidth */
  int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output : super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated left rotations */
  int ldq,          /* leading dimension of Q */
  float *PT,        /* Accumulated right rotations */
  int ldpt,         /* leading dimension of PT */
  float *C,         /* Accumulated left rotations */
  int ldc,          /* leading dimension of C */
  float *work,      /* workspace : ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;

void piro_band_dsbtrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  char *uplo,       /* uplo[0] = 'U' | 'L' for upper | lower triangular stored*/
  int n,            /* #columns in AB */
  int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output : super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated rotations */
  int ldq,          /* leading dimension of Q */
  double *work,     /* workspace : ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;

void piro_band_ssbtrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  char *uplo,       /* uplo[0] = 'U' | 'L' for upper | lower triangular stored*/
  int n,            /* #columns in AB */
  int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output : super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated rotations */
  int ldq,          /* leading dimension of Q */
  float *work,      /* workspace : ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;

void piro_band_zhbtrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  char *uplo,       /* uplo[0] = 'U' | 'L' for upper | lower triangular stored*/
  int n,            /* #columns in AB */
  int ku,           /* upper bandwidth */
  double *AB,       /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  double *D,        /* output : main diagonal of the bidiagonal matrix */
  double *E,        /* output : super diagonal of the bidiagonal matrix */
  double *Q,        /* Accumulated rotations */
  int ldq,          /* leading dimension of Q */
  double *work,     /* workspace :  ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;

void piro_band_chbtrd
(
  char *vect,       /* vect[0] = 'B' | 'Q' | 'P' for both Q and PT | Q | PT  */
  char *uplo,       /* uplo[0] = 'U' | 'L' for upper | lower triangular stored*/
  int n,            /* #columns in AB */
  int ku,           /* upper bandwidth */
  float *AB,        /* input matrix in packed band format */
  int ldab,         /* leading dimension of AB */
  float *D,         /* output : main diagonal of the bidiagonal matrix */
  float *E,         /* output : super diagonal of the bidiagonal matrix */
  float *Q,         /* Accumulated rotations */
  int ldq,          /* leading dimension of Q */
  float *work,      /* workspace : ignored now */
  int *info         /* info = 0 for success and negative for failure */
) ;


void piro_band_dgbbrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  SuiteSparse_long m,   /* #rows in AB */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ncc, /* #columns in C */
  SuiteSparse_long kl,  /* lower bandwidth */
  SuiteSparse_long ku,  /* upper bandwidth */
  double *AB,           /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  double *D,            /* output : main diagonal of the bidiagonal matrix */
  double *E,            /* output : super diagonal of the bidiagonal matrix */
  double *Q,            /* Accumulated left rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  double *PT,           /* Accumulated right rotations */
  SuiteSparse_long ldpt,/* leading dimension of PT */
  double *C,            /* Accumulated left rotations */
  SuiteSparse_long ldc, /* leading dimension of C */
  double *work,         /* workspace : ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

void piro_band_sgbbrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  SuiteSparse_long m,   /* #rows in AB */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ncc, /* #columns in C */
  SuiteSparse_long kl,  /* lower bandwidth */
  SuiteSparse_long ku,  /* upper bandwidth */
  float *AB,            /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  float *D,             /* output : main diagonal of the bidiagonal matrix */
  float *E,             /* output : super diagonal of the bidiagonal matrix */
  float *Q,             /* Accumulated left rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  float *PT,            /* Accumulated right rotations */
  SuiteSparse_long ldpt,/* leading dimension of PT */
  float *C,             /* Accumulated left rotations */
  SuiteSparse_long ldc, /* leading dimension of C */
  float *work,          /* workspace : ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

void piro_band_zgbbrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  SuiteSparse_long m,   /* #rows in AB */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ncc, /* #columns in C */
  SuiteSparse_long kl,  /* lower bandwidth */
  SuiteSparse_long ku,  /* upper bandwidth */
  double *AB,           /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  double *D,            /* output : main diagonal of the bidiagonal matrix */
  double *E,            /* output : super diagonal of the bidiagonal matrix */
  double *Q,            /* Accumulated left rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  double *PT,           /* Accumulated right rotations */
  SuiteSparse_long ldpt,/* leading dimension of PT */
  double *C,            /* Accumulated left rotations */
  SuiteSparse_long ldc, /* leading dimension of C */
  double *work,         /* workspace : ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

void piro_band_cgbbrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  SuiteSparse_long m,   /* #rows in AB */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ncc, /* #columns in C */
  SuiteSparse_long kl,  /* lower bandwidth */
  SuiteSparse_long ku,  /* upper bandwidth */
  float *AB,            /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  float *D,             /* output : main diagonal of the bidiagonal matrix */
  float *E,             /* output : super diagonal of the bidiagonal matrix */
  float *Q,             /* Accumulated left rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  float *PT,            /* Accumulated right rotations */
  SuiteSparse_long ldpt,/* leading dimension of PT */
  float *C,             /* Accumulated left rotations */
  SuiteSparse_long ldc, /* leading dimension of C */
  float *work,          /* workspace : ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

void piro_band_dsbtrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  char *uplo,           /* uplo[0] = 'U' | 'L' for upper | lower tri. stored */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ku,  /* upper bandwidth */
  double *AB,           /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  double *D,            /* output : main diagonal of the bidiagonal matrix */
  double *E,            /* output : super diagonal of the bidiagonal matrix */
  double *Q,            /* Accumulated rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  double *work,         /* workspace : ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

void piro_band_ssbtrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  char *uplo,           /* uplo[0] = 'U' | 'L' for upper | lower tri. stored */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ku,  /* upper bandwidth */
  float *AB,            /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  float *D,             /* output : main diagonal of the bidiagonal matrix */
  float *E,             /* output : super diagonal of the bidiagonal matrix */
  float *Q,             /* Accumulated rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  float *work,          /* workspace : ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

void piro_band_zhbtrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  char *uplo,           /* uplo[0] = 'U' | 'L' for upper | lower tri stored */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ku,  /* upper bandwidth */
  double *AB,           /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  double *D,            /* output : main diagonal of the bidiagonal matrix */
  double *E,            /* output : super diagonal of the bidiagonal matrix */
  double *Q,            /* Accumulated rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  double *work,         /* workspace :  ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

void piro_band_chbtrd_l
(
  char *vect,           /* vect[0] = 'B' | 'Q' | 'P'
                           for both Q and PT | Q | PT  */
  char *uplo,           /* uplo[0] = 'U' | 'L' for upper | lower tri. stored */
  SuiteSparse_long n,   /* #columns in AB */
  SuiteSparse_long ku,  /* upper bandwidth */
  float *AB,            /* input matrix in packed band format */
  SuiteSparse_long ldab,/* leading dimension of AB */
  float *D,             /* output : main diagonal of the bidiagonal matrix */
  float *E,             /* output : super diagonal of the bidiagonal matrix */
  float *Q,             /* Accumulated rotations */
  SuiteSparse_long ldq, /* leading dimension of Q */
  float *work,          /* workspace : ignored now */
  SuiteSparse_long *info/* info = 0 for success and negative for failure */
) ;

#ifdef __cplusplus
}
#endif

#endif /* PIRO_BAND_LAPACK_H */
