/* ========================================================================== */
/* === PIRO_BAND/Test/piro_band_lapackmex.c ================================= */
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
 * Mex function for testing the LAPACK style interface for piro_band methods.
 * See Test/piro_band_lapack.m for details.
 */

#include "piro_band_matlab.h"
#include "piro_band_lapack.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    Int m, n, minmn, i, j, bl, bu, nc, nr, ncl, nrl, work, ldu, ldv,
        crow, err, iscomplex, msize, b ;
    Int *Ap, *Ai ;
    double *Ax, *Axi, *Aband, *b1, *b2, *U, *VT, *rU, *rUi, *rV, *rVi ;
    mxArray *Bmat ;
    char vect [1] ;
    piro_band_mx_options opts ;

    if (nargin < 1 || nargout > 3)
    {
        mexErrMsgTxt ("Usage: [B,U,V]=piro_band_lapack(A,opts)") ;
    }

    if (!mxIsSparse (pargin [0]))
    {
        mexErrMsgTxt ("Input matrix must be sparse") ;
    }

    /* get input parameters */
    piro_bandmex_get_opts (pargin, nargin, 1, &opts) ;

    if (opts.sym && nargout > 2)
    {
        mexErrMsgTxt ("Too many output arguments for symmetric case") ;
    }

    n = mxGetN (pargin [0]) ;
    m = mxGetM (pargin [0]) ;
    minmn = MIN (m, n) ;
    iscomplex = (Int) mxIsComplex (pargin [0]) ;

    if (n < 2)
    {
        mexErrMsgTxt ("size(A,1) must be 2 or more") ;
    }

    /* Allocate space for U */
    vect [0] = 'N' ;
    if (nargout > 1)
    {
        msize = piro_bandmex_multiply (m,m) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        U = (double *) mxMalloc (MAX (1,msize)) ;
        ldu = m ;
        vect [0] = opts.sym ? 'V' : 'Q' ;
    }
    else
    {
        U = NULL ;
        ldu = 0 ;
    }

    /* Allocate space for V (only used for the unsymmetric case) */
    if (nargout > 2)
    {
        msize = piro_bandmex_multiply (n,n) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        VT = (double *) mxMalloc (MAX (1,msize)) ;
        ldv = n ;
        vect [0] = 'B' ;
    }
    else
    {
        VT = NULL ;
        ldv = 0 ;
    }

    /* Allocate space for the bidiagonals */
    msize = minmn ;
    b1 = (double *) mxMalloc (MAX (1,msize) * sizeof (double)) ;
    b2 = (double *) mxMalloc (MAX (1,msize) * sizeof (double)) ;

    Ax = mxGetPr (pargin [0]) ;
    Axi = NULL ;
    if (iscomplex)
    {
        Axi = mxGetPi (pargin [0]) ;
    }

    /* Find the bandwidth of the sparse matrix and store it in packed band
     * format */
    Ap = (Int *) mxGetJc (pargin [0]) ;
    Ai = (Int *) mxGetIr (pargin [0]) ;
    piro_bandmex_find_bandwidth (m, n, Ap, Ai, &bl, &bu) ;

    if (opts.sym)
    {
        if (opts.uplo == 'U')
        {
            /* ignore the lower triangular part of A */
            bl = 0 ;
        }
        else
        {
            /* ignore the upper triangular part of A */
            bu = 0 ;
        }
    }

    crow = bl+bu+1 ;

    msize = piro_bandmex_multiply (crow, n) ;
    msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof(double)) ;
    Aband = (double *) mxCalloc (MAX (1,msize), 1) ;

    piro_bandmex_storeband (Ap, Ai, Ax, Axi, m, n, Aband, bu, bl,
        opts.sym, opts.uplo) ;

    /* reduce the band matrix to the bidiagonal form using LAPACK style
     * inteface. */
    if (minmn > 0) b2 [minmn-1] = 0.0 ;
    if (!opts.sym)
    {
        if (iscomplex)
        {
            piro_band_zgbbrd_l (vect, m, n, 0, bl, bu, Aband, crow,
                            b1, b2, U, ldu, VT, ldv, NULL, 0, NULL, &err) ;
        }
        else
        {
            piro_band_dgbbrd_l (vect, m, n, 0, bl, bu, Aband, crow,
                            b1, b2, U, ldu, VT, ldv, NULL, 0, NULL, &err) ;
        }
    }
    else
    {
        b = (opts.uplo == 'L') ? bl : bu ;
        if (iscomplex)
        {
            piro_band_zhbtrd_l (vect, &(opts.uplo), n, b, Aband, crow,
                            b1, b2, U, ldu, NULL, &err) ;
        }
        else
        {
            piro_band_dsbtrd_l (vect, &(opts.uplo), n, b, Aband, crow,
                            b1, b2, U, ldu, NULL, &err) ;
        }
    }

    if (err != 0)
    {
        piro_bandmex_error (err) ;
    }

    /* free workspace */
    mxFree (Aband) ;

    /* return B as a sparse matrix (bidiagonal or tridiagonal) */
    pargout [0] = piro_bandmex_create_bidiagonal (m, n, opts.sym ? 3:2,
        b1, b2, -1) ;
    mxFree (b1) ;
    mxFree (b2) ;

    /* return U to MATLAB */
    if (nargout > 1)
    {
        pargout [1] = piro_bandmex_put_dense (m, m, &U, iscomplex, 0) ;
    }

    /* return V = VT' to MATLAB */
    if (nargout > 2)
    {
        pargout [2] = piro_bandmex_put_dense (n, n, &VT, iscomplex, 1) ;
    }
}
