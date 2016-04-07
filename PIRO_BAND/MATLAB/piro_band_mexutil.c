/* ========================================================================== */
/* === PIRO_BAND/MATLAB/piro_band_mexutil.c ================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* General utility functions for MATLAB routines
 * */

#include "piro_band_matlab.h"
#include "piro_band.h"
#include "piro_band_internal.h"

/* ========================================================================== */
/* === piro_bandmex_storeband =============================================== */
/* ========================================================================== */

/* Store a sparse matrix in packed band format for a given upper and lower
 * bandwidth. We assume that the Cx array is initialized to zero and fill those
 * entries that are nonzero in the sparse band matrix.
 *
 * If sym is true, then only half of the matrix is considered on input:
 * uplo='L' means the lower part is used and the upper part is assumed to be
 * the transpose of the lower part.  uplo='U' is the reverse.
 * */
void piro_bandmex_storeband
(
    Int *Ap,            /* TODO comment me (and all files) */
    Int *Ai,            /* TODO comment me */
    double *Ax,         /* TODO comment me */
    double *Axi,        /* TODO comment me */
    Int m,              /* TODO comment me */
    Int n,              /* TODO comment me */
    double *Cx,         /* TODO comment me */
    Int bu,             /* TODO comment me */
    Int bl,             /* TODO comment me */
    Int sym,            /* TODO comment me */
    char uplo           /* TODO comment me */
)
{
    Int k ;
    Int p ;
    Int i ;
    Int index ;
    Int obu, ldab ;
    Int islower = (uplo == 'L') ;

    obu = bu;
    ldab = bu + bl + 1 ;
    for (k = 0 ; k < n ; k++)
    {
        for (p = Ap[k] ; p < Ap[k+1] ; p++)
        {
            i = Ai [p] ;
            if (sym && (islower ? (i < k) : (i > k)))
            {
                /* skip this entry */
                continue ;
            }
            index = INDEX (i, k) ;
            if (Axi != NULL)
            {
                Cx[2*index] = Ax[p] ;
                Cx[2*index+1] = Axi[p] ;
            }
            else
            {
                Cx[index] = Ax[p] ;
            }
        }
    }
}

/* ========================================================================== */
/* === piro_bandmex_find_bandwidth ========================================== */
/* ========================================================================== */

/* Find the bandwidth of a sparse matrix using the data structure not the
 * numerical values */
void piro_bandmex_find_bandwidth
(
    Int m,
    Int n,
    Int *Ap,
    Int *Ai,
    Int *bl,
    Int *bu
)
{
    Int k ;
    Int tbl, tbu ;

    tbl = 0 ;
    tbu = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        if (Ap[k] == Ap[k+1]) continue ;
        tbu = MAX(k-Ai[Ap[k]], tbu) ;
        tbl = MAX(Ai[Ap[k+1]-1]-k, tbl) ;
    }
    *bl = tbl ;
    *bu = tbu ;
}

/* ========================================================================== */
/* === piro_bandmex_find_full_bandwidth ===================================== */
/* ========================================================================== */

/* Find the bandwidth of a full matrix using the numerical values. Need to call
 * this function twice for complex matrices (once each with real and imaginary
 * values)
 * */
void piro_bandmex_find_full_bandwidth
(
    Int m,
    Int n,
    double *Ax,
    Int *bl,
    Int *bu
)
{
    Int k, i ;
    Int tbl, tbu ;

    tbu = 0 ;
    tbl = 0 ;
    for (k = 0 ; k < MIN(m, n) ; k++)
    {
        for (i = n-1 ; i > k ; i--)
        {
            /* row access */
            if (Ax[i*m+k] != 0.0) break ;
        }
        tbu = MAX(i - k, tbu) ;

        for (i = m-1 ; i > k ; i--)
        {
            if (Ax[k*m+i] != 0.0) break ;
        }
        tbl = MAX(i - k, tbl) ;
    }
    *bl = MAX(*bl, tbl) ;
    *bu = MAX(*bu, tbu) ;
}

/* ========================================================================== */
/* === piro_bandmex_storeband_withzeroes_full =============================== */
/* ========================================================================== */

/* Store a full matrix in packed band format for a given upper and lower
 * bandwidth
 * */
void piro_bandmex_storeband_withzeroes_full
(
    double *Ax,
    double *Axi,
    Int m,
    Int n,
    double *Cx,
    Int ub,
    Int lb
)
{
    Int k ;
    Int p ;
    Int cnnz ;
    Int i, index ;

    if (m == 0 || n == 0) return ;

    cnnz = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        if (k < ub)
        {
            /*
            for ( i = 0 ; i < ub-k ; i++)
            {
                if (Axi != NULL)
                {
                    Cx[2*cnnz] = 0.0 ;
                    Cx[2*cnnz+1] = 0.0 ;
                }
                else
                {
                    Cx[cnnz] = 0.0 ;
                }
                cnnz++ ;
            }
            */
            cnnz += (ub-k) ;
        }

        for (p = MAX((k*m)+k-ub, k*m) ; p < MIN(k*m+k+lb+1, (k+1)*m) ; p++)
        {
            if (Axi != NULL)
            {
                Cx[2*cnnz] = Ax[p] ;
                Cx[2*cnnz+1] = Axi[p] ;
            }
            else
            {
                Cx[cnnz] = Ax[p] ;
            }
            cnnz++ ;
        }

        index = (k*m+k+lb) - ((k+1)*m) ;
        if (k >= m-lb)
        {
            for (  ; index >= 0 ; index--)
            {
            /*
                if (Axi != NULL)
                {
                    Cx[2*cnnz] = 0.0 ;
                    Cx[2*cnnz+1] = 0.0 ;
                }
                else
                {
                    Cx[cnnz] = 0.0 ;
                }
            */
                cnnz++ ;
            }
        }
    }
}

/* ========================================================================== */
/* === piro_bandmex_band_conjugate_transpose ================================ */
/* ========================================================================== */

/* Find the conjugate transpose of a given band matrix */
void piro_bandmex_band_conjugate_transpose
(
    Int m,
    Int n,
    Int bl,
    Int bu,
    double *A,
    double *AT,
    int iscomplex
)
{
    Int c1, c2 ;
    Int ldab ;
    Int i, j ;
    Int sindex, eindex ;

    ldab = bl + bu ; /*actually it is bl + bu + 1, we will adjust below */

    for (j = 0 ; j < n ; j++)
    {
        c1 = ldab * j + bu ;    /* part of INDEX, -j missing because of ldab */
        c2 = bl + j ;           /* part of INDEX, will skip -i in i loop     */

        sindex = MAX(0, j - bu) ;
        eindex = MIN(m-1, j + bl) ;
        for (i = sindex ; i <= eindex ; i++)
        {
            if (iscomplex)
            {
                AT[2*(c2+ldab*i)] = A[2*(i+c1)] ;       /* rest of INDEX */
                AT[2*(c2+ldab*i)+1] = -A[2*(i+c1)+1] ;  /* rest of INDEX */
            }
            else
            {
                AT[c2+ldab*i] = A[i+c1] ;               /* rest of INDEX */
            }
        }
    }
}

/* ========================================================================== */
/* === piro_bandmex_create_bidiagonal ======================================= */
/* ========================================================================== */

/* create a sparse diagonal, bidiagonal, or symmetric tridiagonal matrix.
   only works for real matrices. */

mxArray *piro_bandmex_create_bidiagonal     /* returns the new sparse matrix */
(
    Int m,          /* number of rows */
    Int n,          /* number of columns */
    Int kind,       /* 1: diagonal, 2: bidiagonal, 3: symmetric tridiagonal */
    double *b1,     /* diagonal entries, length min(m,n) */
    double *b2,     /* offidagonal entries, length min(m,n) */
    Int offset      /* -1: first offdiagonal entry is b2[0]
                        0: first offdiagonal entry is b2[1] */
)
{
    Int mn, bnz, i, j, *Bp, *Bi ;
    mxArray *B ;
    double *Bx ;

    mn = MIN (m,n) ;
    B = mxCreateSparse (m, n, 1 + kind * mn, mxREAL) ;
    Bp = (Int *) mxGetJc (B) ;
    Bi = (Int *) mxGetIr (B) ;
    Bx = mxGetPr (B) ;

    bnz = 0 ;
    for (j = 0 ; j < mn ; j++)
    {
        Bp [j] = bnz ;
        for (i = ((kind > 1)?(j-1):j) ; i <= ((kind < 3) ? j:(j+1)) ; i++)
        {
            if (i >= 0 && i < mn)
            {
                Bi [bnz] = i ;
                if (i < j)
                {
                    Bx [bnz] = b2 [j+offset] ;      /* entry above diagonal */
                }
                else if (i == j)
                {
                    Bx [bnz] = b1 [j] ;             /* diagonal entry */
                }
                else /* i > j */
                {
                    Bx [bnz] = b2 [j+1+offset] ;    /* entry below diagonal */
                }
                if (Bx [bnz] != 0) bnz++ ;          /* drop zeros */
            }
        }
    }
    for ( ; j <= n ; j++)
    {
        Bp [j] = bnz ;
    }
    return (B) ;
}

/* ========================================================================== */
/* === piro_bandmex_multiply ================================================ */
/* ========================================================================== */

/* multiply two integers and check for integer overflow */
Int piro_bandmex_multiply
(
    Int a,
    Int b
)
{
    Int y = a*b ;
    double x = ((double) a) * ((double) b) ;
    if (x != (double) y)
    {
        mexErrMsgTxt ("Problem too large") ;
    }
    return (y) ;
}

/* ========================================================================== */
/* === piro_bandmex_add ===================================================== */
/* ========================================================================== */

/* add two integers and check for integer overflow */
Int piro_bandmex_add
(
    Int a,
    Int b
)
{
    Int y = a+b ;
    double x = ((double) a) + ((double) b) ;
    if (x != (double) y)
    {
        mexErrMsgTxt ("Problem too large") ;
    }
    return (y) ;
}

/* ========================================================================== */
/* === piro_bandmex_error =================================================== */
/* ========================================================================== */

/* print an error message */

void piro_bandmex_error (Int err)
{
    switch (err)
    {
        case PIRO_BAND_OK:
            break ;

        case PIRO_BAND_BLKSIZE_INVALID:
            mexErrMsgTxt ("blocksize invalid") ;
            break ;

        case PIRO_BAND_BL_INVALID:
            mexErrMsgTxt ("lower bandwidth invalid") ;
            break ;

        case PIRO_BAND_BU_INVALID:
            mexErrMsgTxt ("upper bandwidth invalid") ;
            break ;

        case PIRO_BAND_OUT_OF_MEMORY:
            mexErrMsgTxt ("out of memory") ;
            break ;

        case PIRO_BAND_LAPACK_FAILURE:
            mexWarnMsgIdAndTxt ("MATLAB:svd:svdNoConvergence",
                "PIRO_BAND SVD did not converge.") ;
            break ;

        default:
            mexPrintf ("error code: %g\n", (double) err) ;
            mexErrMsgTxt ("internal error!") ;
            break ;
    }
}

/* ========================================================================== */
/* === piro_bandmex_put_dense =============================================== */
/* ========================================================================== */

/* Return a dense matrix to MATLAB, optionally transposing the result.
   A shallow copy is made if X is real and not transposed, or if real,
   transposed, and square.  Otherwise, a copy is made and then X is freed.
 */

mxArray *piro_bandmex_put_dense     /* returns the new MATLAB array */
(
    Int m,
    Int n,
    double **Xhandle,   /* X is m-by-n */
    int iscomplex,      /* true if X is complex, false if X is real */
    int transpose       /* true if X' is to be returned */
)
{
    mxArray *A ;
    Int i, j ;
    double *Areal, *Aimag, *X ;
    X = *Xhandle ;
    if (X == NULL) mexErrMsgTxt ("X is null!") ;

    if (iscomplex)
    {

        /* ------------------------------------------------------------------ */
        /* complex case */
        /* ------------------------------------------------------------------ */

        if (transpose)
        {
            /* return the complex conjugate transpose */
            A = mxCreateDoubleMatrix (n, m, mxCOMPLEX) ;
            Areal = mxGetPr (A) ;
            Aimag = mxGetPi (A) ;
            for (j = 0 ; j < n ; j++)
            {
                for (i = 0 ; i < m ; i++)
                {
                    Areal [j+i*n] =  X [2*(i+j*m)] ;
                    Aimag [j+i*n] = -X [2*(i+j*m)+1] ;
                }
            }
        }
        else
        {
            /* return X unmodified, but split into real / imaginary parts */
            A = mxCreateDoubleMatrix (m, n, mxCOMPLEX) ;
            Areal = mxGetPr (A) ;
            Aimag = mxGetPi (A) ;
            for (j = 0 ; j < n ; j++)
            {
                for (i = 0 ; i < m ; i++)
                {
                    Areal [i+j*m] = X [2*(i+j*m)] ;
                    Aimag [i+j*m] = X [2*(i+j*m)+1] ;
                }
            }
        }
        mxFree (X) ;

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* real case */
        /* ------------------------------------------------------------------ */

        A = NULL ;
        if (transpose)
        {
            if (m == n)
            {
                /* do an in-place transpose, then a shallow copy (below) */
                piro_band_inplace_conjugate_transpose_drl (m, X, m) ;
            }
            else
            {
                /* do an out-of-place transpose and then free X */
                A = mxCreateDoubleMatrix (n, m, mxREAL) ;
                Areal = mxGetPr (A) ;
                for (j = 0 ; j < n ; j++)
                {
                    for (i = 0 ; i < m ; i++)
                    {
                        Areal [j+i*n] =  X [i+j*m] ;
                    }
                }
                mxFree (X) ;
            }
        }
        if (A == NULL)
        {
            /* make a shallow copy */
            A = mxCreateDoubleMatrix (0, 0, mxREAL) ;
            mxSetM (A, m) ;
            mxSetN (A, n) ;
            if (mxGetPr (A) != NULL) mxFree (mxGetPr (A)) ;
            mxSetPr (A, X) ;
        }
    }

    /* set X to NULL to tell the caller not to free it.  It has now already
       been freed, or it is now part of the mxArray returned to MATLAB.  */
    (*Xhandle) = NULL ;

    return (A) ;
}


/* ========================================================================== */
/* === piro_bandmex_getarg ================================================== */
/* ========================================================================== */

/* get a single scalar (0 or 1) from the opts struct */

static Int piro_bandmex_getarg
(
    const mxArray *mxopts,  /* the MATLAB options struct */
    char *field             /* which field to extract from mxopts */
)
{
    int f ;
    mxArray *arg ;

    f = mxGetFieldNumber (mxopts, field) ;
    /* printf ("field: %s %d\n", field, f) ; */
    if (f >= 0)
    {
        arg = mxGetFieldByNumber (mxopts, 0, f) ;
        if (!mxIsEmpty (arg) && (mxIsNumeric (arg) || mxIsLogical (arg)))
        {
            return ((Int) mxGetScalar (arg)) ;
        }
        else
        {
            mexPrintf ("opts.%s:\n", field) ;
            mexErrMsgTxt ("option must be a non-empty numeric value") ;
        }
    }
    return (0) ;        /* the default, if the field is not present */
}


/* ========================================================================== */
/* === piro_bandmex_get_blocksizes ========================================== */
/* ========================================================================== */

/* get the blocksizes from an array of size 4 */

static void piro_bandmex_get_blocksizes
(
    const mxArray *arg,         /* mxArray with the block sizes */
    piro_band_mx_options *opts
)
{
    double *x ;
    int i, j ;

    x = mxGetPr (arg) ;
    for (i = 0 ; i < 4 ; i++)
    {
        if (x [i] < 0)
        {
            /* if any block size is < 0, use defaults for all */
            for (j = 0 ; j < 4 ; j++)
            {
                opts->blks [j] = -1 ;
            }
            return ;
        }
        opts->blks [i] = (Int) x [i] ;
    }
}


/* ========================================================================== */
/* === piro_bandmex_get_opts ================================================ */
/* ========================================================================== */

void piro_bandmex_get_opts
(
    /* inputs */
    const mxArray *pargin [ ],  /* pargin from mexFunction */
    int nargin,                 /* nargin form mexFunction */
    int oparg,                  /* pargin [oparg] is the first option arg */

    /* output */
    piro_band_mx_options *opts
)
{
    Int i, j, k, f ;
    char *s ;
    const mxArray *mxopts, *arg ;
    double *x ;

    /* ---------------------------------------------------------------------- */
    /* set the default options */
    /* ---------------------------------------------------------------------- */

    opts->sym = 0 ;         /* assume matrix is unsymmetric */
    opts->uplo = 'U' ;      /* if symmetric, use the upper triangular part */
    opts->econ = 0 ;        /* full SVD, not economy */
    opts->benchmark = 0 ;   /* do not compute the flop count */
    for (i = 0 ; i < 4 ; i++)
    {
        opts->blks [i] = -1 ;   /* if < 0, default block size is found later */
    }

    if (oparg >= nargin)
    {
        return ;            /* no options passed in; use defaults */
    }

    /* ---------------------------------------------------------------------- */
    /* get the options */
    /* ---------------------------------------------------------------------- */

    if (mxIsStruct (pargin [oparg]))
    {

        /* ------------------------------------------------------------------ */
        /* extract the options from the opts struct */
        /* ------------------------------------------------------------------ */

        /*
            opts.sym        0 (unsymmetric) or non-zero (symmetric)
            opts.uplo       'upper' or 'lower'
            opts.econ       0 (full SVD) or 1 (economy SVD)
            opts.benchmark  0 (no flop counts) or 1 (compute flop counts)
            opts.blks       a double vector of length 4
        */

        if (nargin > oparg + 1)
        {
            mexErrMsgTxt ("opts struct must be final argument") ;
        }
        mxopts = pargin [oparg] ;

        opts->sym  = piro_bandmex_getarg (mxopts, "sym") != 0 ;
        opts->econ = piro_bandmex_getarg (mxopts, "econ") != 0 ;
        opts->benchmark = piro_bandmex_getarg (mxopts, "benchmark") != 0 ;

        /* ------------------------------------------------------------------ */
        /* get opts.uplo */
        /* ------------------------------------------------------------------ */

        f = mxGetFieldNumber (mxopts, "uplo") ;
        if (f >= 0)
        {
            arg = mxGetFieldByNumber (mxopts, 0, f) ;
            if (mxIsChar (arg))
            {
                s = mxArrayToString (arg) ;
                if (strcmp (s, "lower") == 0)
                {
                    opts->uplo = 'L' ;
                }
                else if (strcmp (s, "upper") == 0)
                {
                    opts->uplo = 'U' ;
                }
                else
                {
                    mexPrintf ("opts.uplo = '%s'\n", s) ;
                    mexErrMsgTxt ("unrecognized opts.uplo option") ;
                }
                mxFree (s) ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* get opts.blks */
        /* ------------------------------------------------------------------ */

        f = mxGetFieldNumber (mxopts, "blks") ;
        if (f >= 0)
        {
            arg = mxGetFieldByNumber (mxopts, 0, f) ;
            if (mxGetNumberOfElements (arg) == 4 && mxIsDouble (arg))
            {
                piro_bandmex_get_blocksizes (arg, opts) ;
            }
        }

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* look for options passed as individual arguments */
        /* ------------------------------------------------------------------ */

        /*
            'sym':      set opts.sym = 1
            'upper':    set opts.uplo = 'U' (the default)
            'lower':    set opts.uplo = 'L'
            an array of size 4: use as the opts->blks array

            Note that opts.benchmark cannot be passed in with this method.
        */

        for (k = oparg ; k < nargin ; k++)
        {
            arg = pargin [k] ;
            if (mxIsChar (arg))
            {
                s = mxArrayToString (arg) ;
                if (strcmp (s, "sym") == 0)
                {
                    opts->sym = 1 ;
                }
                else if (strcmp (s, "lower") == 0)
                {
                    opts->uplo = 'L' ;
                }
                else if (strcmp (s, "upper") == 0)
                {
                    opts->uplo = 'U' ;      /* the default */
                }
                else if (strcmp (s, "econ") == 0)
                {
                    opts->econ = 1 ;        /* the economy SVD */
                }
                else
                {
                    mexPrintf ("argument %d: '%s'\n", k+1, s) ;
                    mexErrMsgTxt ("unrecognized option") ;
                }
                mxFree (s) ;
            }
            else if (mxGetNumberOfElements (arg) == 4 && mxIsDouble (arg))
            {
                piro_bandmex_get_blocksizes (arg, opts) ;
            }
            else
            {
                mexPrintf ("argument %d:\n", k+1) ;
                mexErrMsgTxt ("unrecognized option") ;
            }
        }
    }
}

/* ========================================================================== */
/* === piro_bandmex_blocksizes ============================================== */
/* ========================================================================== */

void piro_bandmex_blocksizes
(
    Int rc,
    Int cc,
    Int bl,
    Int bu,
    Int wantuv,

    /* input/output */
    piro_band_mx_options *opts
)
{
    if (opts->blks [0] < 0)
    {
        /* Find the block size for the bidiagonal reduction */
        PIRO_BAND_LONG_NAME (get_blocksize) (rc, cc, bl, bu, wantuv,
            opts->blks) ;
    }

    /* make sure blocksizes are sane */

    if (bl == 0)
    {
        /* no workspace needed for lower triangular part */
        opts->blks [2] = 0 ;
        opts->blks [3] = 0 ;
    }
    else
    {
        opts->blks [2] = MAX (opts->blks [2], (bl > 1 ? 1:0)) ;
        opts->blks [3] = MAX (opts->blks [3], (bl > 1 ? 1:0)) ;
        opts->blks [2] = MIN (opts->blks [2], bl) ;
        opts->blks [3] = MIN (opts->blks [3], bl) ;
        if (opts->blks [2] + opts->blks [3] >= bl)
        {
            opts->blks [2] = bl ;
            opts->blks [3] = 1 ;
        }
    }

    if (bu == 0)
    {
        /* no workspace needed for upper triangular part */
        opts->blks [0] = 0 ;
        opts->blks [1] = 0 ;
    }
    else
    {
        opts->blks [0] = MAX (opts->blks [0], (bu > 1 ? 1:0)) ;
        opts->blks [1] = MAX (opts->blks [1], (bu > 1 ? 1:0)) ;
        opts->blks [0] = MIN (opts->blks [0], bu-1) ;
        opts->blks [1] = MIN (opts->blks [1], bu-1) ;
        if (opts->blks [0] + opts->blks [1] >= bu)
        {
            opts->blks [0] = bu - 1 ;
            opts->blks [1] = 1 ;
        }
    }
}


/* ========================================================================== */
/* === piro_bandmex_identity ================================================ */
/* ========================================================================== */

/* X = eye (n) */

mxArray *piro_bandmex_identity
(
    Int n
)
{
    mxArray *X ;
    double *x ;
    Int k ;
    X = mxCreateDoubleMatrix (n, n, mxREAL) ;
    x = mxGetPr (X) ;
    for (k = 0 ; k < n ; k++)
    {
        x [k*(n+1)] = 1.0 ;
    }
    return (X) ;
}

/* ========================================================================== */
/* === piro_bandmex_scalar ================================================== */
/* ========================================================================== */

/* x = xr + 1i*xi */

mxArray *piro_bandmex_scalar
(
    double xr,      /* real part */
    double xi       /* imaginary part */
)
{
    mxArray *X ;
    double *x ;
    if (xi == 0)
    {
        X = mxCreateDoubleScalar (xr) ;
    }
    else
    {
        X = mxCreateDoubleMatrix (1, 1, mxCOMPLEX) ;
        x = mxGetPr (X) ;
        x [0] = xr ;
        x = mxGetPi (X) ;
        x [0] = xi ;
    }
    return (X) ;
}
