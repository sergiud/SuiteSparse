/* ========================================================================== */
/* === PIRO_BAND/MATLAB/qr_unsqueezemex.c =================================== */
/* ========================================================================== */

/* Returns a row permutation vector that unsqueeze's the 'squeezed R' from
 * the MATLAB sparse QR.  R must be upper triangular or squeezed upper
 * triangular and m-by-n with m >= n.  R must also be sparse.
 *
 * Copyright 2012, Timothy A. Davis, University of Florida (TODO ... !)
 */

#include "mex.h"
#define Int mwSignedIndex

void mexFunction
(
    int nargout,
    mxArray *pargout [],
    int nargin,
    const mxArray *pargin []
)
{
    Int m, n, i, k, rank, *Rp, *Ri ;
    double *perm ;

    /* get inputs */
    if (nargin != 1 || nargout > 2)
    {
        mexErrMsgTxt ("Usage: [p, r] = qr_unsqueeze (R)") ;
    }
    if (!mxIsSparse (pargin [0]))
    {
        mexErrMsgTxt ("R must be sparse") ;
    }

    m = mxGetM (pargin [0]) ;
    n = mxGetN (pargin [0]) ;
    Rp = (Int *) mxGetJc (pargin [0]) ;
    Ri = (Int *) mxGetIr (pargin [0]) ;

    if (m < n)
    {
        mexErrMsgTxt ("R must be m-by-n with m >= n") ;
    }

    /* create perm output */
    pargout [0] = mxCreateDoubleMatrix (1, m, mxREAL) ;
    perm = mxGetPr (pargout [0]) ;

    /* look for live columns in R */
    rank = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        /* find the row index of the last entry in column k, if any */
        i = (Rp [k+1] > Rp [k]) ?  (Ri [Rp [k+1] - 1]) : (-1) ;
        if (i > rank)
        {
            mexErrMsgTxt ("R is not a squeezed upper triangular matrix") ;
        }
        else if (i == rank)
        {
            /* the "diagonal" exists, as entry R (i,k) */
            perm [k] = rank++ ;
        }
        else
        {
            /* column k is dead.  Column k has no diagonal entry */
            perm [k] = -1 ;
        }
    }

    /* create second output, the rank of R as estimated by qr(A) */
    if (nargout > 1)
    {
        pargout [1] = mxCreateDoubleScalar ((double) rank) ;
    }

    /* fill-in any holes in the permutation */
    for (k = 0 ; k < n ; k++)
    {
        if (perm [k] == -1)
        {
            perm [k] = rank++ ;
        }
    }
    for ( ; k < m ; k++)
    {
        perm [k] = k ;
    }

    /* shift by one for MATLAB */
    for (k = 0 ; k < m ; k++)
    {
        perm [k]++ ;
    }
}
