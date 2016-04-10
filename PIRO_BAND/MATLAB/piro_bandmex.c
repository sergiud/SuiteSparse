/* ========================================================================== */
/* === PIRO_BAND/MATLAB/piro_bandmex.c ====================================== */
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
 * Usage :
 * [B, U, V] = piro_band (A)
 * or
 * [B, U, V] = piro_band (A, opts)
 *
 * See piro_band.m for details.
 */

#include "SuiteSparse_config.h"
#include "piro_band.h"
#include "piro_band_matlab.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [],
    int nargin,
    const mxArray *pargin []
)
{
    mxArray *Bmat ;
    Int *Ap, *Ai, *Bp, *Bi ;
    double *Ax, *Axi, *Bx, *dws, *Aband, *b1, *b2, *U, *V ;
    Int m, n, mn, i, j, bnz, bl, bu, nc, nr, ncl, nrl, work, ldu, ldv, crow,
        work1, work2, err, iscomplex, msize ;
    piro_band_mx_options opts ;

    if (nargin < 1 || nargout > 3)
    {
        mexErrMsgTxt ("Usage: [B,U,V] = piro_band (A,options)") ;
    }

    /* get input parameters */
    piro_bandmex_get_opts (pargin, nargin, 1, &opts) ;
    if (opts.uplo != 'U')
    {
        mexErrMsgTxt ("opts.uplo must be 'upper'") ;
    }

    if (opts.benchmark)
    {
        /* flop count will be returned, so the effective nargout is
           decremented, so that U and V are returned properly */
        nargout-- ;
    }

    n = mxGetN (pargin [0]) ;
    m = mxGetM (pargin [0]) ;
    mn = MIN (m, n) ;
    iscomplex = (Int) mxIsComplex (pargin [0]) ;

    if (n < 2)
    {
        mexErrMsgTxt ("size(A,1) must be 2 or more") ;
    }

    if (m != n && opts.sym)
    {
        mexErrMsgTxt ("Invalid input: Matrix can't be rectangular & symmetric");
    }

    /* Allocate space for U.  Freed, or kept, by piro_bandmex_put_dense */
    if (nargout > 1)
    {
        msize = piro_bandmex_multiply (m,m) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        U = (double *) mxMalloc (MAX (1, msize)) ;
        ldu = m ;
    }
    else
    {
        U = NULL ;
        ldu = 0 ;
    }

    /* Allocate space for U.  Freed, or kept, by piro_bandmex_put_dense */
    if (nargout > 2)
    {
        msize = piro_bandmex_multiply (n,n) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        V = (double *) mxMalloc (MAX (1, msize)) ;
        ldv = n ;
    }
    else
    {
        V = NULL ;
        ldv = 0 ;
    }

    /* Allocate space for the bidiagonals */
    b1 = (double *) mxMalloc (MAX (1,mn) * sizeof (double)) ;
    b2 = (double *) mxMalloc (MAX (1,mn) * sizeof (double)) ;

    Ax = mxGetPr (pargin [0]) ;
    Axi = NULL ;
    if (iscomplex)
    {
        Axi = mxGetPi (pargin [0]) ;
    }

    if (mxIsSparse (pargin [0]))
    {
        /* Find the bandwidth of the sparse matrix and store it in packed band
         * format */
        Ap = (Int *) mxGetJc (pargin [0]) ;
        Ai = (Int *) mxGetIr (pargin [0]) ;
        piro_bandmex_find_bandwidth (m, n, Ap, Ai, &bl, &bu) ;
    }
    else
    {
        /* Find the bandwidth of the full matrix and store it in packed band
         * format */
        bl = 0 ;
        bu = 0 ;
        piro_bandmex_find_full_bandwidth (m, n, Ax, &bl, &bu) ;
        if (iscomplex)
        {
            piro_bandmex_find_full_bandwidth (m, n, Axi, &bl, &bu) ;
        }
    }

    /* Force bu to be at least 1 */
    if (bu == 0)
    {
        bu = 1 ;
    }

    if (opts.sym)
    {
        /* ignore the lower triangular part of A */
        bl = 0 ;
    }

    crow = bl+bu+1 ;
    msize = piro_bandmex_multiply (crow,n) ;
    msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
    /* Initialize the band here to simplifiy handling of zeros within the
     * band */
    Aband = (double *) mxCalloc (MAX (1,msize), 1) ;

    if (mxIsSparse (pargin [0]))
    {
        piro_bandmex_storeband (Ap, Ai, Ax, Axi, m, n, Aband, bu, bl,
            opts.sym, 'U') ;
    }
    else
    {
        piro_bandmex_storeband_withzeroes_full (Ax, Axi, m, n, Aband, bu, bl) ;
    }

    piro_bandmex_blocksizes (m, n, bl, bu, (U != NULL || V != NULL), &opts) ;

    /* Allocate workspace */
    work1 = piro_bandmex_multiply (opts.blks [0], opts.blks [1]) ;
    work2 = piro_bandmex_multiply (opts.blks [2], opts.blks [3]) ;
    work = MAX (work1, work2) ;

    /* 2 double values for each column and row rotation */
    msize = piro_bandmex_multiply ((iscomplex ? 4:2) * sizeof (double), work) ;
    dws = (double *) mxMalloc (MAX (1,msize)) ;

    /* reduce the band matrix to the bidiagonal form */
    PIRO_BAND_CLEAR_FLOPS ;
    b2 [0] = 0 ;
    if (iscomplex)
    {
        err = piro_band_reduce_dcl (opts.blks, m, n, 0, bl, bu, Aband, crow,
                        b1, b2+1, U, ldu, V, ldv, NULL, 0, dws, opts.sym) ;
    }
    else
    {
        err = piro_band_reduce_drl (opts.blks, m, n, 0, bl, bu, Aband, crow,
                        b1, b2+1, U, ldu, V, ldv, NULL, 0, dws, opts.sym) ;
    }

    /* free workspace */
    mxFree (dws) ;
    mxFree (Aband) ;

    if (err != 0)
    {
        piro_bandmex_error (err) ;
    }

    /* return U to MATLAB */
    if (nargout > 1)
    {
        pargout [1] = piro_bandmex_put_dense (m, m, &U, iscomplex, 0) ;
    }

    /* return V to MATLAB */
    if (nargout > 2)
    {
        pargout [2] = piro_bandmex_put_dense (n, n, &V, iscomplex, 0) ;
    }

    /* return B as a sparse matrix (bidiagonal or tridiagonal) */
    pargout [0] = piro_bandmex_create_bidiagonal (m, n, opts.sym ? 3:2,
        b1, b2, 0) ;

    /* return the flop count as the last argument */
    if (opts.benchmark)
    {
        pargout [nargout] = mxCreateDoubleScalar (PIRO_BAND_GET_FLOPS) ;
    }

    /* Free workspace */
    mxFree (b1) ;
    mxFree (b2) ;
}
