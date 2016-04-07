/* ========================================================================== */
/* === PIRO_BAND/Test/piro_band_qrmex.c ===================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Mex function for QR factorization of a matrix. Uses the left looking
 * band QR.   See piro_band_qr.m for details.  No input options are used.
 */

#include "piro_band_matlab.h"
#include "piro_band.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [],
    int nargin,
    const mxArray *pargin []
)
{
    Int m, n, bl, bu, crow, i, j, k, err, ldq, ldr, ldv, dsize, bandi, nentries,
        computeQ, nnz, msize, minmn, ri1, ri2, obu, ldab, rindex, iscomplex,
        issparse ;
    Int *Ap, *Ai, *Rp, *Ri ;
    double *Ax, *Axi, *dws, *Aband, *Q, *R, *V, *X1, *work, *beta, *tmp,
        *rU, *rUi, *Rreal, *Rimag, *rbeta, *rbetai ;
    mxArray *mV, *mBeta, *mR ;
    static const char *QRnames [ ] = { "V", "beta", "R", "bl", "bu", "m", "n" };

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (nargin != 1 || nargout > 2)
    {
        mexErrMsgTxt ("Usage: [Q,R]=piro_band_qr (A) or QR=piro_band_qr (A)") ;
    }

    computeQ = (nargout == 2) ;

    /* ---------------------------------------------------------------------- */
    /* get input matrix A */
    /* ---------------------------------------------------------------------- */

    n = mxGetN (pargin [0]) ;
    m = mxGetM (pargin [0]) ;
    minmn = MIN (m ,n) ;

    iscomplex = (Int) mxIsComplex (pargin [0]) ;
    issparse = (Int) mxIsSparse (pargin [0]) ;
    Ax = mxGetPr (pargin [0]) ;
    Axi = NULL ;
    if (iscomplex)
    {
        Axi = mxGetPi (pargin [0]) ;
    }

    /* ---------------------------------------------------------------------- */
    /* initialize the output matrix Q */
    /* ---------------------------------------------------------------------- */

    if (computeQ)
    {
        msize = piro_bandmex_multiply (m, minmn) ;
        msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
        Q = (double *) mxCalloc (MAX (1,msize), 1) ;
        ldq = m ;
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
        ldq = 0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* Find the bandwidth and store in packed band format */
    /* ---------------------------------------------------------------------- */

    if (issparse)
    {
        Ap = (Int *) mxGetJc (pargin [0]) ;
        Ai = (Int *) mxGetIr (pargin [0]) ;
        piro_bandmex_find_bandwidth (m, n, Ap, Ai, &bl, &bu) ;
    }
    else
    {
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

    crow = bl + bu + 1 ;
    msize = piro_bandmex_multiply (crow, n) ;
    msize = piro_bandmex_multiply (msize, (iscomplex?2:1)*sizeof (double)) ;
    Aband = (double *) mxCalloc (MAX (1,msize), 1) ;

    if (issparse)
    {
        /* the sparse matrix A is not symmetric; both upper/lower parts used */
        piro_bandmex_storeband (Ap, Ai, Ax, Axi, m, n, Aband, bu, bl, 0, 'X') ;
    }
    else
    {
        piro_bandmex_storeband_withzeroes_full (Ax, Axi, m, n, Aband, bu, bl) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Allocate workspace for the QR factorization */
    /* ---------------------------------------------------------------------- */

    /* 2 * bl + bu + 1 workspace for a column   */
    /* (bl + 1) * n for storing the V matrix */
    /* n for the beta */
    dsize = piro_bandmex_multiply (bl, n+2) ;
    dsize = piro_bandmex_add (dsize, bu+1) ;
    dsize = piro_bandmex_add (dsize, 2*n) ;

    if (computeQ)
    {
        /* n for X1 (to compute Q) */
        dsize = piro_bandmex_add (dsize, n) ;
    }

    msize = piro_bandmex_multiply (dsize, (iscomplex ? 2:1) * sizeof (double)) ;
    dws = (double *) mxMalloc (MAX (1,msize)) ;
    tmp = dws ;

    work = tmp ;
    msize = (iscomplex ? 2:1) * ((2 * bl) + bu + 1) ;
    tmp += msize ;

    beta = tmp ;
    msize = (iscomplex ? 2:1) * n ;
    tmp += msize ;

    V = tmp ;
    msize = (iscomplex ? 2:1) * ((bl+1) * n) ;
    tmp += msize ;

    X1 = tmp ; /* not allocated if !computeQ */

    ldv = bl + 1 ;

    /* ---------------------------------------------------------------------- */
    /* Compute the factorization */
    /* ---------------------------------------------------------------------- */

    if (iscomplex)
    {
        err = piro_band_qr_dcl (m, n, bl, bu, Aband, crow, V, ldv, beta, work) ;
    }
    else
    {
        err = piro_band_qr_drl (m, n, bl, bu, Aband, crow, V, ldv, beta, work) ;
    }
    if (err != 0)
    {
        piro_bandmex_error (err) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Compute Q from the householder vectors V in C itself */
    /* ---------------------------------------------------------------------- */

    if (computeQ)
    {
        if (iscomplex)
        {
            piro_band_computeQ_dcl (m, n, bl, V, ldv, beta, m, minmn, Q, ldq,
                            work, X1) ;
        }
        else
        {
            piro_band_computeQ_drl (m, n, bl, V, ldv, beta, m, minmn, Q, ldq,
                            work, X1) ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* return results to MATLAB */
    /* ---------------------------------------------------------------------- */

    if (!computeQ)
    {

        /* ------------------------------------------------------------------ */
        /* QR = piro_band_qr (A) ; return the QR struct */
        /* ------------------------------------------------------------------ */

        if (iscomplex)
        {
            mV = mxCreateDoubleMatrix (bl+1, n, mxCOMPLEX) ;
            mBeta = mxCreateDoubleMatrix (1, n, mxCOMPLEX) ;
            mR = mxCreateDoubleMatrix (bl+bu+1, n, mxCOMPLEX) ;
            rUi = mxGetPi (mV) ;
            rbetai = mxGetPi (mBeta) ;
            Rimag = mxGetPi (mR) ;
        }
        else
        {
            mV = mxCreateDoubleMatrix (bl+1, n, mxREAL) ;
            mBeta = mxCreateDoubleMatrix (1, n, mxREAL) ;
            mR = mxCreateDoubleMatrix (bl+bu+1, n, mxREAL) ;
            rUi = NULL ;
            rbetai = NULL ;
            Rimag = NULL ;
        }

        rU = mxGetPr (mV) ;
        rbeta = mxGetPr (mBeta) ;
        Rreal = mxGetPr (mR) ;

        /* copy V to rU and rUi */
        for (j = 0 ; j < n ; j++)
        {
            for (i = 0 ; i < bl+1 ; i++)
            {
                if (iscomplex)
                {
                    rU  [i+j*(bl+1)] = V [(i+j*(bl+1))*2] ;
                    rUi [i+j*(bl+1)] = V [(i+j*(bl+1))*2+1] ;
                }
                else
                {
                    rU [i+j*(bl+1)] = V [i+j*(bl+1)] ;
                }
            }
        }

        /* copy beta to rbeta and rbetai */
        for (j = 0 ; j < n ; j++)
        {
            if (iscomplex)
            {
                rbeta  [j] = beta [j*2] ;
                rbetai [j] = beta [j*2+1] ;
            }
            else
            {
                rbeta [j] = beta [j] ;
            }
        }

        /* copy R to Rreal and Rimag */
        for (j = 0 ; j < n ; j++)
        {
            for (i = 0 ; i < crow ; i++)
            {
                if (iscomplex)
                {
                    Rreal [i+j*crow] = Aband [(i+j*crow)*2] ;
                    Rimag [i+j*crow] = Aband [(i+j*crow)*2+1] ;
                }
                else
                {
                    Rreal [i+j*crow] = Aband [i+j*crow] ;
                }
            }
        }

        /* Create the structure to return and set all the fields */
        pargout [0] = mxCreateStructMatrix (1, 1, 7, QRnames) ;
        mxSetFieldByNumber (pargout [0], 0, 0, mV) ;
        mxSetFieldByNumber (pargout [0], 0, 1, mBeta) ;
        mxSetFieldByNumber (pargout [0], 0, 2, mR) ;
        mxSetFieldByNumber (pargout [0], 0, 3, mxCreateDoubleScalar (bl)) ;
        mxSetFieldByNumber (pargout [0], 0, 4, mxCreateDoubleScalar (bu)) ;
        mxSetFieldByNumber (pargout [0], 0, 5, mxCreateDoubleScalar (m)) ;
        mxSetFieldByNumber (pargout [0], 0, 6, mxCreateDoubleScalar (n)) ;

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* [Q,R] = piro_band_qr (A) */
        /* ------------------------------------------------------------------ */

        /* return Q to MATLAB */
        pargout [0] = piro_bandmex_put_dense (m, minmn, &Q, iscomplex, 0) ;

        /* return R to MATLAB */

        obu = bl + bu ;
        ldab = obu + 1 ;
        if (issparse)
        {
            if (iscomplex)
            {
                pargout [1] = mxCreateSparse (minmn, n, crow*n, mxCOMPLEX) ;
                Rimag = mxGetPi (pargout [1]) ;
            }
            else
            {
                pargout [1] = mxCreateSparse (minmn, n, crow*n, mxREAL) ;
                Rimag = NULL ;
            }
            Rp = (Int *) mxGetJc (pargout [1]) ;
            Ri = (Int *) mxGetIr (pargout [1]) ;
        }
        else
        {
            if (iscomplex)
            {
                pargout [1] = mxCreateDoubleMatrix (minmn, n, mxCOMPLEX) ;
                Rimag = mxGetPi (pargout [1]) ;
            }
            else
            {
                pargout [1] = mxCreateDoubleMatrix (minmn, n, mxREAL) ;
                Rimag = NULL ;
            }
            Rp = NULL ;
            Ri = NULL ;
        }
        Rreal = mxGetPr (pargout [1]) ;

        nnz = 0 ;
        for (j = 0 ; j < MIN (m+bu, n) ; j++)
        {
            ri1 = MAX ((j - obu), 0) ; /* ri1 <= m-1 */
            ri2 = MIN (j, m-1) ;
            if (issparse)
            {
                Rp [j] = nnz ;
            }

            for (i = ri1 ; i <= ri2 ; i++)
            {
                rindex = i - j + obu + ldab * j ;
                if (iscomplex)
                {
                    if (issparse)
                    {
                        Ri [nnz] = i ;
                        Rreal [nnz] = Aband [rindex*2] ;
                        Rimag [nnz] = Aband [rindex*2+1] ;
                    }
                    else
                    {
                        Rreal [j*m+i] = Aband [rindex*2] ;
                        Rimag [j*m+i] = Aband [rindex*2+1] ;
                    }
                }
                else
                {
                    if (issparse)
                    {
                        Ri [nnz] = i ;
                        Rreal [nnz] = Aband [rindex] ;
                    }
                    else
                    {
                        Rreal [j*m+i] = Aband [rindex] ;
                    }
                }
                /* drop zeros from sparse R */
                if (issparse)
                {
                    if (Rreal [nnz] != 0 || (iscomplex && Rimag [nnz] != 0))
                    {
                        nnz++ ;
                    }
                }
            }
        }
        if (issparse)
        {
            for ( ; j <= n ; j++)
            {
                Rp [j] = nnz ;
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* Free work space */
    /* ---------------------------------------------------------------------- */

    if (Q != NULL) mxFree (Q) ;
    mxFree (Aband) ;
    mxFree (dws) ;
}
