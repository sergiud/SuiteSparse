/* User-includable file */

/* TODO needs LOTS of comments */

#ifndef BLKSKY_H
#define BLKSKY_H

typedef struct sky_common_struct
{
    SuiteSparse_long opt ;
    SuiteSparse_long cspace ;
    SuiteSparse_long elbow ;
    SuiteSparse_long rcount ;

    /* Includes reallocs for A. Not for the work space which may not be there
     * in the final version.
     * */
    SuiteSparse_long realloc_cnt ;
    SuiteSparse_long *iwork ;
    double *dwork ;
    double safemin2 ;
    double safemx2 ;
    /* Estimated bandwidth for R */
    SuiteSparse_long bw_est ;
} sky_common ;

typedef struct sky_symbolic_struct
{
    SuiteSparse_long *trow ;      /* TODO: how big is it?  What is it? etc... */
    SuiteSparse_long *encol ;
    SuiteSparse_long *prev ;
    SuiteSparse_long *next ;
    SuiteSparse_long *min_trow ;  /* TODO what is min_trow? */
    SuiteSparse_long n ;
    double giv_cnt ;
    double swap_cnt ;
    double flp_cnt ;
} sky_symbolic ;

typedef struct sky_skyline_struct
{
    SuiteSparse_long m ;
    SuiteSparse_long n ;
    SuiteSparse_long *cp ;
    double *A ;
} sky_skyline ;

typedef struct sky_factor_struct
{
    double *b1 ;
    double *b2 ;
    double *U ;
    double *V ;
} sky_factor ;

SuiteSparse_long sky_set_common
(
    sky_common *scom 
) ;

SuiteSparse_long sky_allocate_symbolic
(
    sky_common *scom,
    sky_symbolic *sym
) ;

SuiteSparse_long sky_free_symbolic
(
    sky_symbolic *sym
) ;

SuiteSparse_long sky_allocate_skyline
(
    sky_common *scom,
    sky_skyline *sky,
    sky_symbolic *sym
) ;

SuiteSparse_long sky_free_skyline
(
    sky_skyline *sky
) ;

void blksky_find_symbolic
(
    sky_common *scom,
    SuiteSparse_long *Ap,
    SuiteSparse_long *Ai,
    double *Ax,
    SuiteSparse_long m,
    SuiteSparse_long n,
    sky_symbolic *sym
) ;

SuiteSparse_long sky_sparse_to_skyline
(
    sky_skyline *sky,
    SuiteSparse_long *Ap,
    SuiteSparse_long *Ai,
    double *Ax
) ;

SuiteSparse_long sky_reduce
(
    sky_common *scom,
    sky_symbolic *sym,
    sky_skyline *sky,
    sky_factor *sfact
) ;

void givens                 /* TODO rename */
(
    double *a,              /* Entry array of size 1/2, i/p to rotate */
    SuiteSparse_long ai,
    double *b,	            /* Entry array of size 1/2, i/p to rotate */
    SuiteSparse_long bi,
    double *c_op,           /* cosine of the givens rotation */
    SuiteSparse_long ci,
    double *s_op,           /* sine of the givens rotation */
    SuiteSparse_long si,
    double *r_op,           /* rotated vector's first entry */
    SuiteSparse_long ri,
    sky_common *scom
) ;

void apply_givens           /* TODO rename */
(
    double *da,
    double *db,
    double c,
    double s
) ;

#endif
