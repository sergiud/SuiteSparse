/*
 * Reduces sparse upper triangular A, into bidiagonal matrix B, using 
 * rightmost block, leftmost corner method.
 *
 * Usage :
 *
 * B = blksky_mex(A, opt, bw) 
 *
 * opt = 1 => Fixed block size, vanilla symbolic factorization and a static 
 *            skyline. This is the default option if opt is not present.
 * opt = 2 => Dynamic block size, vanilla symbolic factorization and a static 
 *            skyline.
 * opt = 3 => Fixed block size, no symbolic factorization and a dynamically
 *            allocated skyline.
 */

#include "mex.h"
#include "genband_util.h"
#include "genband_internal.h"
#include "blksky.h"
#undef ASSERT
/*
#include "piro_band_internal.h"
*/
#include "piro_band_blas.h"
#include "piro_band_matlab.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    Int m, n, csize, nnz, i, j ;
    Int *Ap, *Ai, *Bp, *Bi ;
    double *Ax, *Bx, *rb1, *rb2 ;
    sky_common sc, *scom ;
    sky_symbolic sm, *ssym ;
    sky_skyline ss, *ssky ;
    sky_factor sf, *sfact ;

    if (nargout > 3 || nargin > 3 || nargin < 1)
    {
	mexErrMsgTxt("Usage: [B,U,V] = blksky (A, opt, bw)") ;
    }

    m = (Int) mxGetM (pargin [0]) ;
    n = (Int) mxGetN (pargin [0]) ;
    /* printf ("m "ID" n "ID" sizeof(Int) "ID"\n", m, n, sizeof (Int)) ; */

    scom = &sc ;
    sky_set_common(scom) ;

    ssym = &sm ;
    ssym->n = n ;

    ssky = &ss ;
    ssky->m = m ;
    ssky->n = n ;

    /*
    printf ("allocating sfact b1 and b2:\n") ;
    printf ("m %g n %g sizeof(double) %d MIN(m,n) %g\n",
        (double) m,
        (double) n,
        sizeof (double),
        (double) MIN(m,n)) ;
    */

    sfact = &sf ;
    sfact->b1 = (double *) mxMalloc(MIN(m, n) * sizeof(double)) ;
    sfact->b2 = (double *) mxMalloc(MIN(m, n) * sizeof(double)) ;

    /*
    printf ("sfact %p %p\n", sfact->b1, sfact->b2) ;
    */

    if (nargout > 1)
    {
        /* TODO document this */
        pargout [1] = mxCreateDoubleMatrix(m, m, mxREAL) ;
        sfact->U = mxGetPr(pargout [1]) ;
        for (i = 0 ; i < m ; i++)
        {
            sfact->U[i*m+i] = 1.0 ;
        }
    }
    else
    {
        sfact->U = NULL ;
    }

    if (nargout > 2)
    {
        /* TODO document this */
        pargout [2] = mxCreateDoubleMatrix(n, n, mxREAL) ;
        sfact->V = mxGetPr(pargout [2]) ;
        for (i = 0 ; i < n ; i++)
        {
            sfact->V[i*n+i] = 1.0 ;
        }
    }
    else
    {
        sfact->V = NULL ;
    }

    Ap = (Int *) mxGetJc (pargin [0]) ;
    Ai = (Int *) mxGetIr (pargin [0]) ;
    Ax = mxGetPr (pargin [0]) ;
    nnz = Ap [n] ;

    /* TODO use a struct here */
    if (nargin < 2)
    {
        scom->opt = 1 ;
    }
    else
    {
        scom->opt = (Int) mxGetScalar (pargin [1]) ;
    }

    if (nargin < 3)
    {
        scom->bw_est = 1 ;
    }
    else
    {
        double x ;
        x = mxGetScalar (pargin [2]) ;
        scom->bw_est = (Int) x ;
    }

    scom->bw_est = MAX (1, scom->bw_est) ;
    scom->bw_est = MIN (n, scom->bw_est) ;

    if (scom->opt == 3)
    {
        scom->cspace = 5 ;
        scom->elbow = 3 * n ;
    }
    else
    {
        scom->cspace = 0 ;
        scom->elbow = 0 ;
    }

    sky_allocate_symbolic (scom, ssym) ;
    blksky_find_symbolic (scom, Ap, Ai, Ax, m, n, ssym) ;
    sky_allocate_skyline (scom, ssky, ssym) ;
    sky_sparse_to_skyline (ssky, Ap, Ai, Ax) ;

    if (scom->opt != 2)
    {
        scom->rcount = 32 ;
        scom->dwork = (double *) mxMalloc( 4 * scom->rcount * sizeof(double)) ;
        scom->iwork = (Int *) mxMalloc( 2 * scom->rcount * sizeof(Int)) ;
    }
    else
    {
        scom->rcount = -1 ;
        scom->dwork = NULL ;
        scom->iwork = NULL ;
    }

    sky_reduce (scom, ssym, ssky, sfact) ;

    /* return B as a sparse bidiagonal matrix */
    pargout [0] = piro_bandmex_create_bidiagonal (m, n, 2,
        sfact->b1, sfact->b2, 0) ;

    mxFree(sfact->b1) ;
    mxFree(sfact->b2) ;
    if (scom->opt != 2)
    {
        mxFree(scom->dwork) ;
        mxFree(scom->iwork) ;
    }

    sky_free_skyline(ssky) ; 
    sky_free_symbolic(ssym) ;
}
