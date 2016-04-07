/* ========================================================================== */
/* === PIRO_BAND/MATLAB/piro_band_qr.c ====================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Computes the left looking QR factorization of a band matrix stored in the
 * packed format. The Q is in terms of the householder vectors stored in V and
 * the band matrix A is overwritten with the R.
 */

#include "piro_band_blas.h"
#include "piro_band.h"


/* Apply the householder reflections to work(i1 : i2). The size of v should be
 * appropriate too. */
static void PIRO_BAND (happly)
(
    Entry *v,       /* TODO comment me */
    Int j,
    Int vsize,
    Entry *beta,
    Entry *work,
    Int i1
)
{
    BLAS_INT m, n, ldc, incv ;

    char LR[1] ;
    Entry hwork[CSIZE(1)] ;
    Entry c_beta[CSIZE(1)] ;

    LR[0] = 'L' ;

    /* BLAS integers: */
    m = (BLAS_INT) vsize ;
    n = 1 ;
    incv = 1 ;
    ldc = (BLAS_INT) vsize ;

#ifdef PIROBAND_COMPLEX
    CONJ(c_beta, beta) ;
    LAPACK_ZLARF(LR, &m, &n, v, &incv, c_beta, work+i1, &ldc, hwork) ;
#else
    LAPACK_DLARF(LR, &m, &n, v, &incv, beta, work+i1, &ldc, hwork) ;
#endif

}

/* Generate the householder transformation to reduce the vector x */
static void PIRO_BAND (house)
(
    Entry *x,
    Entry *beta,
    Int n,
    Entry *diag
)
{
    Entry alpha[CSIZE(1)] ;
    BLAS_INT n2, incx ;

    n2 = (BLAS_INT) n ;
    incx = 1 ;

    ASSIGN_TO_SCALAR(alpha, x, 0) ;

#ifdef PIROBAND_COMPLEX
    PRINT_QRVALUES("alpha[%d] =", 0, alpha, 0) ;
    LAPACK_ZLARFG(&n2, alpha, x+2, &incx, beta) ;
    PRINT_QRVALUES("beta[%d] =", 0,  alpha, 0) ;
#else
    PRINT_QRVALUES("alpha[%d] =", 0, alpha, 0) ;
    LAPACK_DLARFG(&n2, alpha, x+1, &incx, beta) ;
    PRINT_QRVALUES("beta[%d] =", 0, alpha, 0) ;
#endif

    ASSIGN_TO_SCALAR(diag, alpha, 0) ;
    x[0] = 1.0 ;
#ifdef PIROBAND_COMPLEX
    x[1] = 0.0 ;
#endif
    return ;
}

#ifndef QR_INDEX
#define QR_INDEX(row, col) ( (row)-(col)+(obu)+(bl))
#endif

/* Find the QR factorization for a band matrix using householder
 * transformations. */
int PIRO_BAND(qr)           /* returns error code */
(
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    Entry *A,               /* Band Matrix                             */
    Int ldab,               /* leading dimension of Ax                 */
    Entry *V,               /* o/p householder vectors                 */
    Int ldv,                /* leading dimension of V                  */
    Entry *Beta,            /* o/p accumulated right rotations         */
    Entry *work             /* workspace of size ((2 * bl) + bu + 1)   */
)
{
    Int i, j, k, width, obu, nz, vsize, i1, i2, info, rindex ;
    Entry beta[CSIZE(1)] ;
    Entry diag[CSIZE(1)] ;
    Entry one[2] ;
    Entry *v ;

    /* Inputs are checked by caller (piro_band_svdkernel):
     *      Requres m > 0, n > 0, bl >= 0, bu > 0, ldab >= bl+bu+1,
     *      ldv >= 0; A, V, work, and Beta non-null.
     */

    nz = 2 * bl + bu + 1 ;
    width = bl + bu ;
    one[0] = 1.0 ;
    one[1] = 0.0 ;

    /* obu and ldab are used by INDEX macro */
    obu = bu ;

    /* Initialize the top of the workspace */
    for (i = 0 ; i < bl ; i++)
    {
        ASSIGN_ZERO_TO_MATRIX(work, CSIZE(i)) ;
    }

    /* Compute the QR one column at a time */
    for (k = 0 ; k < MIN(m+bu, n) ; k++)
    {
        PRINT(("k = %d \n", k)) ;
        for (i = 0 ; i < bl + (MAX(0, k - bu) -  (k - bu)) ; i++)
        {
            ASSIGN_ZERO_TO_MATRIX(work, CSIZE(i)) ;
            PRINT_QRVALUES("work[%d] = ", CSIZE(i), work, CSIZE(i)) ;
        }

        /* Copy the current column into the bottom of the workspace */
        for (rindex = MAX( 0, k - bu) ; rindex <= MIN( k + bl , m - 1) ;
                                                        i++, rindex++)
        {
            ASSIGN_MATRIX_TO_MATRIX(work, CSIZE(i), A, INDEX(rindex, k));
            PRINT_QRVALUES("work[%d] = ", CSIZE(i), work, CSIZE(i)) ;
        }

        /* Apply the pending householders from the previous columns */
        for (j = MAX(k - width , 0) ; j < MIN(m, k) ; j++)
        {
            /* Apply the rotation from the jth column */
            /*printf("Applying rotation from the %d th column \n", j) ;*/
            vsize = MIN(j + bl , m - 1 ) - j + 1 ;
            v = V + CSIZE(j * ldv ) ;
            ASSIGN_TO_SCALAR(beta, Beta, CSIZE(j)) ;

            /* We add bl to the index because of the bl additional space on the
             * top of the workspace */
            i1 = QR_INDEX (j, k) ;
            /*i2 = QR_INDEX (MIN(j+bl, m-1), k) ;*/
            /* Apply the rotation to work (i1 : i2) */
            PIRO_BAND(happly) (v, j, vsize, beta, work, CSIZE(i1)) ;
        }

        vsize = MIN(k + bl , m - 1 ) - k + 1 ;
        i1 = QR_INDEX (MIN(k, m-1) , k ) ;
        i2 = QR_INDEX (MIN(k + bl , m - 1) , k ) ;
        PRINT(("i1 = %d i2 = %d \n", i1, i2)) ;
        for (i = i1 ; i <= i2 ; i++)
        {
            PRINT_QRVALUES("work[%d] = ", CSIZE(i), work, CSIZE(i)) ;
        }
        /* Find the householder vector to zero work (i1 : i2) */
        if (vsize > 0)
        {
            PIRO_BAND (house) (work+CSIZE(i1), Beta + CSIZE(k), vsize, diag)  ;
        }
        else
        {
            ASSIGN_ZERO_TO_MATRIX(Beta, CSIZE(k)) ;
            /*printf("Skipping beta assignment *********\n") ;*/
            /* Assign one to V */
            ASSIGN_TO_MATRIX(one, V, CSIZE(k*ldv)) ;
        }

        /* Copy the result to R (which is A here) and V */
        for (i = 0 ; i <= i1 - 1 ; i++)
        {
            ASSIGN_MATRIX_TO_MATRIX(A, CSIZE(k* ldab + i), work, CSIZE(i));
            PRINT_QRVALUES("Ax[%d] = ",  CSIZE(k * ldab + i), work, CSIZE(i)) ;
        }

        if (m <= n && k > m-1)
        {
            ASSIGN_MATRIX_TO_MATRIX(A, CSIZE(k* ldab + i), work, CSIZE(i));
            i++ ; /* To skip the copy in V when vsize < 1. */
        }
        else
        {
            ASSIGN_TO_MATRIX(diag, A, CSIZE(k* ldab + i));
            PRINT_QRVALUES("Ax[%d] = ",  CSIZE(k * ldab + i), diag, 0) ;
        }

        for (rindex = 0 ; i <= i2 ; i++, rindex++)
        {
            ASSIGN_MATRIX_TO_MATRIX(V, CSIZE(k* ldv + rindex), work, CSIZE(i));
            PRINT_QRVALUES("V[%d] = ", CSIZE(k * ldv + rindex), work, CSIZE(i));
        }
        PRINT_QRVALUES("beta[%d] = ", k, Beta, CSIZE(k)) ;
    }

    return (PIRO_BAND_OK) ;

}

/* Compute the Q of the QR factorization from the housholder vectors V and the
 * Beta.
 * V is a band matrix with only the lower band. V will be of the same size as
 * the orignal band matrix.
 * No permutation is used.
 * */
void PIRO_BAND(computeQ)
(
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int bl,                 /* lower bandwidth                         */
    Entry *V,               /* householder vectors                     */
    Int ldv,                /* leading dimension of V                  */
    Entry *Beta,            /* accumulated right rotations             */
    Int mq,                 /* #rows in Q                              */
    Int nq,                 /* #columns in Q                           */
    Entry *Q,               /* o/p Q of the QR factorization           */
    Int ldq,                /* leading dimension of Q                  */
    Entry *work,            /* workspace of size xxxxxxx               */
    Entry *X1               /* workspace of size xxxxxxx               */
)
{
    Int k ;
    Int qk ;
    Int nentries ;
    Int i, vi, qi ;
    Int vsize ;
    Entry beta[CSIZE(1)] ;
    Entry c_V[CSIZE(1)] ;

    vsize = bl + 1 ;
    for (k = n-1 ; k >= 0 ; k-- )
    {
        if (k + bl < m)
        {
            nentries = vsize ;
        }
        else
        {
            nentries = m - k ;
        }

        if (nentries <= 0) continue ;
        PRINT(("K = %d ************\n", k)) ;

        ASSIGN_TO_SCALAR(beta, Beta, CSIZE(k)) ;
        for (vi = k*ldv, i = 0 ; vi < k*ldv+nentries ; vi++, i++)
        {
            ASSIGN_TO_SCALAR(c_V, V, CSIZE(vi)) ;
            CONJ(c_V, c_V) ;
            MULT_I(work, CSIZE(i), beta, 0, c_V, 0) ;
            PRINT_QRVALUES("c_V[%d] =", 0, c_V, 0) ;
            PRINT_QRVALUES("beta[%d] =", 0, beta, 0) ;
            PRINT_QRVALUES("work[%d] =", CSIZE(i), work, CSIZE(i)) ;
        }

        for (qk = 0 ; qk < nq ; qk++)
        {
            PRINT(("QK1 = %d +++++++++  \n", qk)) ;
            /* Find X1[qk] with a dot product */
            ASSIGN_ZERO_TO_MATRIX(X1, CSIZE(qk)) ;
            for (qi = qk * ldq + k, i = 0 ; qi <= qk * ldq + k + nentries - 1 ;
                                                     qi++, i++)
            {
                MULT_ADD_I(X1, CSIZE(qk), Q, CSIZE(qi), work, CSIZE(i)) ;
                PRINT_QRVALUES("X1[%d] =", CSIZE(qk), X1, CSIZE(qk)) ;
            }
        }

        for (qk = 0 ; qk < nq ; qk++)
        {
            PRINT(("QK2 = %d ---------  \n", qk)) ;
            for (qi = qk * ldq + k, vi = k*ldv ;
                    qi <= qk * ldq + k + nentries - 1 ; qi++, vi++)
            {
                MULT_SUB_I(Q, CSIZE(qi), X1, CSIZE(qk), V, CSIZE(vi)) ;
                PRINT_QRVALUES("Q[%d] =", CSIZE(qi), Q, CSIZE(qi)) ;
            }
            PRINT(("\n")) ;
        }
    }
}
