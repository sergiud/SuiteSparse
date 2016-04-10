/*
 * Symbolic factorization for the method that reduces sparse upper triangular A,
 * into bidiagonal matrix B, using rightmost block, leftmost corner method.
 *
 * Usage : 
 * toprow = blksky_symbolic_mex(R) ;
 * [toprow, no_givens] = bllsky_symbolic_mex(R) ;
 * [toprow, no_givens, no_swaps] = bllsky_symbolic_mex(R) ;
 * [toprow, no_givens, no_swaps, flop_count] = bllsky_symbolic_mex(R) ;
 */

#include "mex.h"
#include "genband_util.h"
#include "genband_internal.h"
#include "blksky.h"

void mexFunction
(
    int nlhs,
    mxArray *plhs[],
    int nrhs,
    const mxArray *prhs[]
)
{
    Int m, n ;                  /* dimensions of A */
    Int i ;                     /* index */
    Int *Ap, *Ai ;
    double *Ax ;
    double *trow ;
    Int *min_trow ;
    sky_common sc, *scom ;
    sky_symbolic sm, *ssym ;    /* TODO sometimes called "sym". Be consistent */

    if (nlhs > 4 || nrhs != 1 )
    {
	mexErrMsgTxt("Invalid number of args to blksky_mex\n") ;
    }

    scom = &sc ;
    scom->opt = 1 ;
    scom->bw_est = 1 ;

    m = (Int) mxGetM(prhs[0]) ;
    n = (Int) mxGetN(prhs[0]) ;
    ssym = &sm ;
    ssym->n = n ;

    Ap = (Int *) mxGetJc(prhs[0]) ;
    Ai = (Int *) mxGetIr(prhs[0]) ;
    Ax = mxGetPr(prhs[0]) ;

    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL) ;
    trow = mxGetPr(plhs[0]) ;

    /*min_trow = (int *) mxMalloc( n * sizeof(int)) ;*/
    sky_allocate_symbolic(scom, ssym) ;

    blksky_find_symbolic(scom, Ap, Ai, Ax, m, n, ssym) ;

    /* TODO what is going on here?  How do min_trow and trow differ? */
    min_trow = ssym->min_trow ;

    /* TODO: Oh my!  What is going on here?  Why is min_trow being copied 
       into trow? */
    for (i = 0 ; i < n ; i++)
    {
        trow[i] = min_trow[i] + 1 ;
    }

    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleScalar(ssym->giv_cnt) ;
    }

    if (nlhs > 2)
    {
        plhs[2] = mxCreateDoubleScalar(ssym->swap_cnt) ;
    }

    if (nlhs > 3)
    {
        plhs[3] = mxCreateDoubleScalar(ssym->flp_cnt) ;
    }

    sky_free_symbolic(ssym) ;
    /*mxFree(min_trow) ;*/

}

