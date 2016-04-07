#include "gp_mex.hpp"

namespace SuiteSparse_Mongoose
{

/* get a MATLAB dense column vector */
double *gp_mex_get_double (Int n, const mxArray *X)
{
    return (mxGetPr (X)) ;
}

/* return a Double vector to MATLAB */
double *gp_mex_put_double (Int n, const double *b, mxArray **X)
{
    double *x ;
    Int k ;
    *X = mxCreateDoubleMatrix (n, 1, mxREAL) ;      /* create x */
    x = mxGetPr (*X) ;
    for (k = 0 ; k < n ; k++) x [k] = b [k] ;       /* copy x = b */
    return (x) ;
}

/* get a MATLAB flint array and convert to Int */
Int *gp_mex_get_int
(
    Int n,
    const mxArray *Imatlab,
    Int *imax,
    Int lo
)
{
    double *p ;
    Int i, k, *C = (Int*) malloc(n * sizeof (Int));

    p = mxGetPr (Imatlab) ;
    *imax = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        i = (Int) p[k];
        C [k] = i - 1 ;
        if (i < lo) mexErrMsgTxt ("index out of bounds") ;
        *imax = MONGOOSE_MAX2 (*imax, i) ;
    }
    return (C) ;
}

/* return an Int array to MATLAB as a flint row vector */
mxArray *gp_mex_put_int(Int *p, Int n, Int offset, Int do_free)
{
    mxArray *X = mxCreateDoubleMatrix (1, n, mxREAL) ;
    double *x = mxGetPr (X) ;
    Int k ;
    for (k = 0 ; k < n ; k++) x [k] = (p ? p [k] : k) + offset ;
    if (do_free) free(p);
    return (X) ;
}

/* return an Int array to MATLAB as a flint row vector */
mxArray *gp_mex_put_logical(bool *p, Int n)
{
    mxArray *X = mxCreateDoubleMatrix(1, n, mxREAL);
    double *x = mxGetPr(X);
    for(Int k = 0; k < n; k++) x[k] = p[k] ? 1.0 : 0.0;
    return X;
}

}