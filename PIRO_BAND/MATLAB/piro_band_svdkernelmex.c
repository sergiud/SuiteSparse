/* ========================================================================== */
/* === PIRO_BAND/MATLAB/piro_band_svdkernelmex.c ============================ */
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
 * Mex function to find the SVD of a band matrix in sparse/full format.
 * See piro_band_svd.m for details.  This function is not meant to be
 * called by the end-user.
 *
 * Note that opts.uplo and opts.benchmark are silently ignored.
 */

#include "piro_band_matlab.h"
#include "piro_band.h"
#include "piro_band_blas.h"

#define PIROBAND_LONG
#define PIROBAND_COMPLEX
#include "piro_band_internal.h"
#undef PIROBAND_LONG
#undef PIROBAND_COMPLEX

void mexFunction
(
    int nargout,
    mxArray *pargout [],
    int nargin,
    const mxArray *pargin []
)
{
    Int m, n, i, j, bl, bu, ldu, ldv, ldv1, ldq, crow, err, minmn, nrq, ncq,
        msize, computeQ, iscomplex, itmp, Vm, Vn, Um, Un, rc, cc, work,
        work1, work2, transposeA ;
    BLAS_INT info, ldc, ldu2, ldvt, n2, ncc, ncvt, nru ;
    Int *Ap, *Ai ;
    double *Ax, *Az, *dws, *Aband, *b1, *b2, *U, *VT, *Q, *Q1, *Mat, *BTx ;
    double **Vhandle, **Uhandle ;
    char uplo [1] ;
    piro_band_mx_options opts ;

    dws = NULL ;
    Aband = NULL ;
    b1 = NULL ;
    b2 = NULL ;
    BTx = NULL ;
    U = NULL ;
    VT = NULL ;
    Q = NULL ;
    Q1 = NULL ;

    if (nargin < 1 || nargout > 3 )
    {
        mexErrMsgTxt ("Usage: [U,S,V]=piro_band_svdmex (A,opts)") ;
    }

    /* get input parameters */
    piro_bandmex_get_opts (pargin, nargin, 1, &opts) ;

    n = mxGetN (pargin [0]) ;
    m = mxGetM (pargin [0]) ;
    minmn = MIN (m, n) ;
    iscomplex = (Int) mxIsComplex (pargin [0]) ;

    if (minmn < 2)
    {
        mexErrMsgTxt ("min(size(A)) must be 2 or more") ;
    }

    Ax = mxGetPr (pargin [0]) ;
    Az = NULL ;
    if (iscomplex)
    {
        Az = mxGetPi (pargin [0]) ;
    }

    /* ---------------------------------------------------------------------- */

    computeQ = (opts.econ && m != n) ;
    transposeA = computeQ && (m < n) ;

    /* ---------------------------------------------------------------------- */

    if (transposeA)
    {
        rc = n ;
        cc = m ;
    }
    else
    {
        rc = m ;
        cc = n ;
    }

    /* Allocate space for U */
    if (nargout > 1 && !computeQ)
    {
        msize = piro_bandmex_multiply (m, m) ;
        msize = MAX (1, msize) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        U = (double *) mxMalloc (msize) ;
        ldu = m ;
        nru = (BLAS_INT) m ;
    }
    else
    {
        U = NULL ;
        ldu = 0 ;
        nru = 0 ;
    }

    /* Allocate space for V */
    if (nargout > 2 || (transposeA && nargout > 1))
    {
        msize = piro_bandmex_multiply (cc, cc) ;
        msize = MAX (1, msize) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        VT = (double *) mxMalloc (msize) ;
        ldv = cc ;
        ncvt = (BLAS_INT) cc ;
    }
    else
    {
        VT = NULL ;
        ldv = 0 ;
        ncvt = 0 ;
    }

    /* Create the matrix Q for the QR factorization. */
    if (computeQ)
    {
        msize = piro_bandmex_multiply (rc, minmn) ;
        msize = MAX (1, msize) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        Q  = (double *) mxCalloc (msize, 1) ;
        Q1 = (double *) mxMalloc (msize) ;
        ldq = rc ;
        nrq = rc ;
        ncq = minmn ;
        for (i = 0 ; i < minmn ; i++)
        {
            if (iscomplex)
            {
                Q [2*(i*ldq+i)] = 1.0 ;
            }
            else
            {
                Q [i*ldq+i] = 1.0 ;
            }
        }
    }
    else
    {
        Q = NULL ;
        Q1 = NULL ;
        ldq = 0 ;
        nrq = 0 ;
        ncq = 0 ;
    }

    if (mxIsSparse (pargin [0]))
    {
        /* Find the bandwidth of the sparse A */
        Ap = (Int *) mxGetJc (pargin [0]) ;
        Ai = (Int *) mxGetIr (pargin [0]) ;
        piro_bandmex_find_bandwidth (m, n, Ap, Ai, &bl, &bu) ;
    }
    else
    {
        /* Find the bandwidth of the dense A */
        bl = 0 ;
        bu = 0 ;
        piro_bandmex_find_full_bandwidth (m, n, Ax, &bl, &bu) ;
        if (iscomplex)
        {
            piro_bandmex_find_full_bandwidth (m, n, Az, &bl, &bu) ;
        }
    }

    /* make sure the upper bandwidth is valid (bu cannot be zero) */
    if (transposeA)
    {
        /* the matrix will be transposed first, so make sure bl > 0 holds */
        if (bl == 0)
        {
            bl = 1 ;
        }
    }
    else
    {
        if (bu == 0)
        {
            bu = 1 ;
        }
    }

    /* store in packed band format */
    crow = bl+bu+1 ;
    msize = piro_bandmex_multiply (crow, n) ;
    msize = MAX (1, msize) ;
    msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
    Aband = (double *) mxCalloc (msize, 1) ;

    if (mxIsSparse (pargin [0]))
    {
        /* the sparse matrix A is not symmetric; both upper/lower parts used */
        piro_bandmex_storeband (Ap, Ai, Ax, Az, m, n, Aband, bu, bl, 0, 'X') ;
    }
    else
    {
        piro_bandmex_storeband_withzeroes_full (Ax, Az, m, n, Aband, bu, bl) ;
    }

    if (transposeA)
    {
        /* transpose the matrix */
        BTx = (double *) mxMalloc (msize) ;
        piro_bandmex_band_conjugate_transpose (m, n, bl, bu, Aband,
            BTx, iscomplex) ;
        Mat = BTx ;
        itmp = bl ;
        bl = bu ;
        bu = itmp ;
    }
    else
    {
        BTx = NULL ;
        Mat = Aband ;
    }

    if (computeQ)
    {
        /* Allocate workspace for the QR factorization */
        /* 2 * bl + bu + 1 workspace for a column   */
        /* (bl + 1) * cc for storing the V1 matrix */
        /* cc for the beta */
        double *Work, *beta, *tmp, *V1, *X1 ;

        msize = piro_bandmex_multiply (bl, cc+2) ;
        msize = piro_bandmex_add (msize, bu+1) ;
        msize = piro_bandmex_add (msize, 2*cc) ;
        if (computeQ)
        {
            /* cc for X1 (to compute Q) */
            msize = piro_bandmex_add (msize, cc) ;
        }
        msize = MAX (1, msize) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        dws = (double *) mxMalloc (msize) ;

        tmp = dws ;

        Work = tmp ;
        msize = (iscomplex ? 2:1) * ((2 * bl) + bu + 1) ;
        tmp += msize ;

        beta = tmp ;
        msize = (iscomplex ? 2:1) * cc ;
        tmp += msize ;

        V1 = tmp ;
        msize = (iscomplex ? 2:1) * ((bl+1) * cc) ;
        tmp += msize ;

        X1 = tmp ;

        ldv1 = bl + 1 ;

        /* Compute the QR factorization */
        if (iscomplex)
        {
            err = piro_band_qr_dcl (rc, cc, bl, bu, Mat, crow, V1, ldv1, beta,
                                    Work) ;
        }
        else
        {
            err = piro_band_qr_drl (rc, cc, bl, bu, Mat, crow, V1, ldv1, beta,
                                    Work) ;
        }
        if (err != 0)
        {
            piro_bandmex_error (err) ;
        }

        /* Compute Q from the householder vectors V1.  */
        if (iscomplex)
        {
            piro_band_computeQ_dcl (rc, cc, bl, V1, ldv1, beta, rc, minmn,
                Q, ldq, Work, X1) ;
        }
        else
        {
            piro_band_computeQ_drl (rc, cc, bl, V1, ldv1, beta, rc, minmn,
                Q, ldq, Work, X1) ;
        }

        /* Free the workspace from QR */
        if (dws != NULL) mxFree (dws) ;
        dws = NULL ;

        /* Adjust bl and bu for upper triangular R */
        /* bu = bl + bu will work correctly, but will not use the efficient
         * blocksizes and will only be an workaround. Need to reassign Mat to
         * do it correctly.
         * */
        if (bl+bu > cc-1)
        {
            msize = iscomplex ? 2 * (bl+bu-(cc-1)) : (bl+bu-(cc-1)) ;
            Mat = Mat + msize ;
            bu = cc-1 ;
        }
        else
        {
            bu = bl + bu ;
        }
        /*bu = bl + bu ; */
        bl = 0 ;

    }

    /* Allocate space for the bidiagonals */
    msize = minmn ;
    msize = MAX (2, msize) ;
    b1 = (double *) mxMalloc (msize * sizeof (double)) ;
    b2 = (double *) mxMalloc (msize * sizeof (double)) ;

    piro_bandmex_blocksizes (rc, cc, bl, bu,
        (U != NULL || VT != NULL || Q != NULL), &opts) ;

    /* Allocate workspace */
    work1 = piro_bandmex_multiply (opts.blks [0], opts.blks [1]) ;
    work2 = piro_bandmex_multiply (opts.blks [2], opts.blks [3]) ;
    work = MAX (work1, work2) ;
    work = MAX (2, work) ;


    /* 2 double values for each column and row rotation */
    msize = piro_bandmex_multiply ((iscomplex ? 4:2) * sizeof (double), work) ;
    dws = (double *) mxMalloc (msize) ;

    /* Reduce to bidiagonal matrix */
    if (iscomplex)
    {
        err = piro_band_reduce_dcl (opts.blks, rc, cc, nrq, bl, bu, Mat, crow,
            b1, b2+1, U, ldu, VT, ldv, Q, nrq, dws, 0) ;
    }
    else
    {
        err = piro_band_reduce_drl (opts.blks, rc, cc, nrq, bl, bu, Mat, crow,
            b1, b2+1, U, ldu, VT, ldv, Q, nrq, dws, 0) ;
    }

    if (err != 0)
    {
        piro_bandmex_error (err) ;
    }

    if (VT != NULL)
    {
        /* Need to transpose VT */
        if (iscomplex)
        {
            piro_band_inplace_conjugate_transpose_dcl (cc, VT, ldv) ;
        }
        else
        {
            piro_band_inplace_conjugate_transpose_drl (cc, VT, ldv) ;
        }
    }

    if (computeQ)
    {
        /* Q1 = Q' */
        if (iscomplex)
        {
            piro_band_general_transpose_dcl (rc, minmn, Q, rc, Q1, minmn) ;
        }
        else
        {
            piro_band_general_transpose_drl (rc, minmn, Q, rc, Q1, minmn) ;
        }
    }
    else
    {
        ncq = 1 ; /* for Fortran interface */
    }

    /* Q no longer needed */
    if (Q != NULL) mxFree (Q) ;
    Q = NULL ;

    if (U == NULL) ldu = 1 ;
    if (VT == NULL) ldv = 1 ;
    if (dws != NULL) mxFree (dws) ;
    dws = NULL ;

    /* TBD : if nargout > 1: uses more space than reqd because of lapack. */
    msize = (iscomplex ? 2:1) * (2 * MAX (1, 4 * minmn)) ;
    msize = piro_bandmex_multiply (msize, sizeof (double)) ;
    dws = (double *) mxMalloc (msize) ;

    n2 = (BLAS_INT) minmn ;
    ncc = (BLAS_INT) nrq ;
    ldvt = (BLAS_INT) ldv ;
    ldu2 = (BLAS_INT) ldu ;
    ldc = (BLAS_INT) ncq ;

    uplo [0] = 'U' ;
    if (iscomplex)
    {
        LAPACK_ZBDSQR (uplo, &n2, &ncvt, &nru, &ncc, b1, b2+1, VT, &ldvt, U,
                &ldu2, Q1, &ldc, dws, &info) ;
    }
    else
    {
        LAPACK_DBDSQR (uplo, &n2, &ncvt, &nru, &ncc, b1, b2+1, VT, &ldvt, U,
                &ldu2, Q1, &ldc, dws, &info) ;
    }

    err = (Int) info ;

    if (err != 0)
    {
        /* if positive, then LAPACK did not converge to the solution.
           if negative, then one or more inputs to LAPACK were invalid. */
        piro_bandmex_error ((err > 0) ?
            PIRO_BAND_LAPACK_FAILURE : PIRO_BAND_LAPACK_INVALID) ;
    }

    /* ---------------------------------------------------------------------- */
    /* free workspace no longer needed */
    /* ---------------------------------------------------------------------- */

    if (dws != NULL) mxFree (dws) ;
    dws = NULL ;

    if (b2 != NULL) mxFree (b2) ;
    b2 = NULL ;

    if (Aband != NULL) mxFree (Aband) ;
    Aband = NULL ;

    if (BTx != NULL) mxFree (BTx) ;
    BTx = NULL ;

    /* ---------------------------------------------------------------------- */
    /* return results to MATLAB */
    /* ---------------------------------------------------------------------- */

    if (computeQ)
    {
        if (!transposeA)
        {
            /* return U = Q1' to MATLAB */
            if (nargout > 1)
            {
                pargout [0] = piro_bandmex_put_dense (minmn, rc, &Q1,
                    iscomplex, 1) ;
            }
            /* return V = VT' to MATLAB */
            if (nargout > 2)
            {
                pargout [2] = piro_bandmex_put_dense (cc, cc, &VT,
                    iscomplex, 1) ;
            }
        }
        else
        {
            /* return U = VT' to MATLAB */
            if (nargout > 1)
            {
                pargout [0] = piro_bandmex_put_dense (cc, cc, &VT,
                    iscomplex, 1) ;
            }
            /* return V = Q1' to MATLAB */
            if (nargout > 2)
            {
                pargout [2] = piro_bandmex_put_dense (minmn, rc, &Q1,
                    iscomplex, 1) ;
            }
        }
    }
    else
    {
        /* return U to MATLAB: the U matrix of size rc-by-rc */
        if (nargout > 1)
        {
            pargout [0] = piro_bandmex_put_dense (rc, rc, &U, iscomplex, 0) ;
        }
        /* return V = VT' to MATLAB */
        if (nargout > 2)
        {
            pargout [2] = piro_bandmex_put_dense (cc, cc, &VT,
                iscomplex, 1) ;
        }
    }

    /* return S to MATLAB */
    if (nargout <= 1)
    {
        pargout [0] = piro_bandmex_put_dense (minmn, 1, &b1, 0, 0) ;
    }
    else
    {
        /* return S as a diagonal sparse matrix */
        Int ms = opts.econ ? minmn : m ;
        Int ns = opts.econ ? minmn : n ;
        pargout [1] = piro_bandmex_create_bidiagonal (ms, ns, 1, b1, NULL, 0) ;
    }

    /* Free workspace */
    if (U  != NULL) mxFree (U)  ;
    if (VT != NULL) mxFree (VT) ;
    if (Q1 != NULL) mxFree (Q1) ;
    if (b1 != NULL) mxFree (b1) ;
}
