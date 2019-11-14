/* ========================================================================== */
/* CXSparse/Include/cs.h file */
/* ========================================================================== */

/* This is the CXSparse/Include/cs.h file.  It has the same name (cs.h) as
   the CSparse/Include/cs.h file.  The 'make install' for SuiteSparse installs
   CXSparse, and this file, instead of CSparse.  The two packages have the same
   cs.h include filename, because CXSparse is a superset of CSparse.  Any user
   program that uses CSparse can rely on CXSparse instead, with no change to the
   user code.  The #include "cs.h" line will work for both versions, in user
   code, and the function names and user-visible typedefs from CSparse all
   appear in CXSparse.  For experimenting and changing the package itself, I
   recommend using CSparse since it's simpler and easier to modify.  For
   using the package in production codes, I recommend CXSparse since it has
   more features (support for complex matrices, and both int and long
   versions).
 */

/* ========================================================================== */

#ifndef _CXS_H
#define _CXS_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#include "cs_export.h"

#ifndef NCOMPLEX
#include <complex.h>
#endif /* !defined(NCOMPLEX) */

#if defined(__cplusplus)
extern "C" {
#endif

#ifndef NCOMPLEX
#if defined(CXSPARSE_C99_COMPLEX)
typedef double _Complex cs_complex_t;

#define CS_COMPLEX_MULCC(x, y) ((x) * (y))
#define CS_COMPLEX_MULCR(x, y) ((x) * (y))
#define CS_COMPLEX_MUL(x, y) ((x) * (y))
#define CS_COMPLEX_MAKE_ENTRY(x) (x)
#define CS_COMPLEX_ZERO() CS_COMPLEX_MAKE_ENTRY(0.)
#define CS_COMPLEX_SQRT(x) csqrt(x)
#define CS_COMPLEX_ABS(x) cabs(x)
#define CS_COMPLEX_ONE() CS_COMPLEX_MAKE_ENTRY(1.)
#define CS_COMPLEX_MINUS_ONE() CS_COMPLEX_MAKE_ENTRY(-1.)
#define CS_COMPLEX_ADD(x, y) ((x) + (y))
#define CS_COMPLEX_SUB(x, y) ((x) - (y))
#define CS_COMPLEX_IS_ZERO(x) ((x) == 0)
#define CS_COMPLEX_DIV(x, y) ((x) / (y))
#define CS_COMPLEX_MAKE(re, im) ((re) + (im) * _Complex_I)

#elif defined(CXSPARSE_C99_COMPLEX_SUBSET)
typedef _Dcomplex cs_complex_t;

#define CS_COMPLEX_MUL(x, y) CS_COMPLEX_MULCC(x, y)
#define CS_COMPLEX_MAKE_ENTRY(x) CS_COMPLEX_MAKE(x, 0.0)
#define CS_COMPLEX_ZERO() CS_COMPLEX_MAKE_ENTRY(0.)
#define CS_COMPLEX_SQRT(x) csqrt(x)
#define CS_COMPLEX_ABS(x) cabs(x)
#define CS_COMPLEX_ONE() CS_COMPLEX_MAKE_ENTRY(1.)
#define CS_COMPLEX_MINUS_ONE() CS_COMPLEX_MAKE_ENTRY(-1.)
#define CS_COMPLEX_ADD(x, y) _DCOMPLEX_(creal(x) + creal(y), cimag(x) + cimag(y))
#define CS_COMPLEX_SUB(x, y) _DCOMPLEX_(creal(x) - creal(y), cimag(x) - cimag(y))
#define CS_COMPLEX_IS_ZERO(x) (creal(x) == 0 && cimag(x) == 0)
#define CS_COMPLEX_DIV(x, y) _DCOMPLEX_((x._Val[0] * y._Val[0] + x._Val[1] * y._Val[1]) / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1]),  (x._Val[1] * y._Val[0] - x._Val[0] * y._Val[1]) / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1]))
#define CS_COMPLEX_MAKE(re, im) _DCOMPLEX_(re, im)
#define CS_COMPLEX_MULCC(x, y) _Cmulcc(x, y)
#define CS_COMPLEX_MULCR(x, y) _Cmulcr(x, y)
#else
#error "C99 complex support not specified"
#endif /* defined(CXSPARSE_C99_COMPLEX) */
#endif /* !defined(NCOMPLEX) */

#ifdef CS_COMPLEX
#ifdef NCOMPLEX
#error "C99 complex support undefined"
#endif /* defined(NCOMPLEX) */
#define CS_MUL(x, y) CS_COMPLEX_MUL(x, y)
#define CS_MAKE_ENTRY(x) CS_COMPLEX_MAKE_ENTRY(x)
#define CS_ZERO() CS_COMPLEX_ZERO()
#define CS_SQRT(x) CS_COMPLEX_SQRT(x)
#define CS_ONE() CS_COMPLEX_ONE()
#define CS_MINUS_ONE() CS_COMPLEX_MINUS_ONE()
#define CS_ADD(x, y) CS_COMPLEX_ADD(x, y)
#define CS_SUB(x, y) CS_COMPLEX_SUB(x, y)
#define CS_IS_ZERO(x) CS_COMPLEX_IS_ZERO(x)
#define CS_DIV(x, y) CS_COMPLEX_DIV(x, y)
#else /* !defined(CS_COMPLEX) */
#define CS_MUL(x, y) ((x) * (y))
#define CS_MAKE_ENTRY(x) (x)
#define CS_SQRT(x) sqrt(x)
#define CS_ZERO() 0
#define CS_ONE() 1
#define CS_MINUS_ONE() (-1)
#define CS_ADD(x, y) ((x) + (y))
#define CS_SUB(x, y) ((x) - (y))
#define CS_IS_ZERO(x) ((x) == 0)
#define CS_DIV(x, y) ((x) / (y))
#endif /* defined(CS_COMPLEX) */

#define CS_VER 3                    /* CXSparse Version */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Sept 12, 2017"       /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2016"
#define CXSPARSE

#include "SuiteSparse_config.h"
#define cs_long_t       SuiteSparse_long
#define cs_long_t_id    SuiteSparse_long_id
#define cs_long_t_max   SuiteSparse_long_max

/* -------------------------------------------------------------------------- */
/* double/int version of CXSparse */
/* -------------------------------------------------------------------------- */

/* --- primary CSparse routines and data structures ------------------------- */

typedef struct cs_di_sparse  /* matrix in compressed-column or triplet form */
{
    int nzmax ;     /* maximum number of entries */
    int m ;         /* number of rows */
    int n ;         /* number of columns */
    int *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    int nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs_di ;

CXSPARSE_EXPORT cs_di *cs_di_add (const cs_di *A, const cs_di *B, double alpha, double beta) ;
CXSPARSE_EXPORT int cs_di_cholsol (int order, const cs_di *A, double *b) ;
CXSPARSE_EXPORT int cs_di_dupl (cs_di *A) ;
CXSPARSE_EXPORT int cs_di_entry (cs_di *T, int i, int j, double x) ;
CXSPARSE_EXPORT int cs_di_lusol (int order, const cs_di *A, double *b, double tol) ;
CXSPARSE_EXPORT int cs_di_gaxpy (const cs_di *A, const double *x, double *y) ;
CXSPARSE_EXPORT cs_di *cs_di_multiply (const cs_di *A, const cs_di *B) ;
CXSPARSE_EXPORT int cs_di_qrsol (int order, const cs_di *A, double *b) ;
CXSPARSE_EXPORT cs_di *cs_di_transpose (const cs_di *A, int values) ;
CXSPARSE_EXPORT cs_di *cs_di_compress (const cs_di *T) ;
CXSPARSE_EXPORT double cs_di_norm (const cs_di *A) ;
CXSPARSE_EXPORT int cs_di_print (const cs_di *A, int brief) ;
CXSPARSE_EXPORT cs_di *cs_di_load (FILE *f) ;

/* utilities */
CXSPARSE_EXPORT void *cs_di_calloc (int n, size_t size) ;
CXSPARSE_EXPORT void *cs_di_free (void *p) ;
CXSPARSE_EXPORT void *cs_di_realloc (void *p, int n, size_t size, int *ok) ;
CXSPARSE_EXPORT cs_di *cs_di_spalloc (int m, int n, int nzmax, int values, int t) ;
CXSPARSE_EXPORT cs_di *cs_di_spfree (cs_di *A) ;
CXSPARSE_EXPORT int cs_di_sprealloc (cs_di *A, int nzmax) ;
CXSPARSE_EXPORT void *cs_di_malloc (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

typedef struct cs_di_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    int *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    int *q ;        /* fill-reducing column permutation for LU and QR */
    int *parent ;   /* elimination tree for Cholesky and QR */
    int *cp ;       /* column pointers for Cholesky, row counts for QR */
    int *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    int m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;    /* # entries in U for LU; in R for QR */
} cs_dis ;

typedef struct cs_di_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs_di *L ;      /* L for LU and Cholesky, V for QR */
    cs_di *U ;      /* U for LU, r for QR, not used for Cholesky */
    int *pinv ;     /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} cs_din ;

typedef struct cs_di_dmperm_results    /* cs_di_dmperm or cs_di_scc output */
{
    int *p ;        /* size m, row permutation */
    int *q ;        /* size n, column permutation */
    int *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    int *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    int nb ;        /* # of blocks in fine dmperm decomposition */
    int rr [5] ;    /* coarse row decomposition */
    int cc [5] ;    /* coarse column decomposition */
} cs_did ;

CXSPARSE_EXPORT int *cs_di_amd (int order, const cs_di *A) ;
CXSPARSE_EXPORT cs_din *cs_di_chol (const cs_di *A, const cs_dis *S) ;
CXSPARSE_EXPORT cs_did *cs_di_dmperm (const cs_di *A, int seed) ;
CXSPARSE_EXPORT int cs_di_droptol (cs_di *A, double tol) ;
CXSPARSE_EXPORT int cs_di_dropzeros (cs_di *A) ;
CXSPARSE_EXPORT int cs_di_happly (const cs_di *V, int i, double beta, double *x) ;
CXSPARSE_EXPORT int cs_di_ipvec (const int *p, const double *b, double *x, int n) ;
CXSPARSE_EXPORT int cs_di_lsolve (const cs_di *L, double *x) ;
CXSPARSE_EXPORT int cs_di_ltsolve (const cs_di *L, double *x) ;
CXSPARSE_EXPORT cs_din *cs_di_lu (const cs_di *A, const cs_dis *S, double tol) ;
CXSPARSE_EXPORT cs_di *cs_di_permute (const cs_di *A, const int *pinv, const int *q,
    int values) ;
CXSPARSE_EXPORT int *cs_di_pinv (const int *p, int n) ;
CXSPARSE_EXPORT int cs_di_pvec (const int *p, const double *b, double *x, int n) ;
CXSPARSE_EXPORT cs_din *cs_di_qr (const cs_di *A, const cs_dis *S) ;
CXSPARSE_EXPORT cs_dis *cs_di_schol (int order, const cs_di *A) ;
CXSPARSE_EXPORT cs_dis *cs_di_sqr (int order, const cs_di *A, int qr) ;
CXSPARSE_EXPORT cs_di *cs_di_symperm (const cs_di *A, const int *pinv, int values) ;
CXSPARSE_EXPORT int cs_di_usolve (const cs_di *U, double *x) ;
CXSPARSE_EXPORT int cs_di_utsolve (const cs_di *U, double *x) ;
CXSPARSE_EXPORT int cs_di_updown (cs_di *L, int sigma, const cs_di *C, const int *parent) ;

/* utilities */
CXSPARSE_EXPORT cs_dis *cs_di_sfree (cs_dis *S) ;
CXSPARSE_EXPORT cs_din *cs_di_nfree (cs_din *N) ;
CXSPARSE_EXPORT cs_did *cs_di_dfree (cs_did *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */

CXSPARSE_EXPORT int *cs_di_counts (const cs_di *A, const int *parent, const int *post,
    int ata) ;
CXSPARSE_EXPORT double cs_di_cumsum (int *p, int *c, int n) ;
CXSPARSE_EXPORT int cs_di_dfs (int j, cs_di *G, int top, int *xi, int *pstack,
    const int *pinv) ;
CXSPARSE_EXPORT int *cs_di_etree (const cs_di *A, int ata) ;
CXSPARSE_EXPORT int cs_di_fkeep (cs_di *A, int (*fkeep) (int, int, double, void *),
    void *other) ;
CXSPARSE_EXPORT double cs_di_house (double *x, double *beta, int n) ;
CXSPARSE_EXPORT int *cs_di_maxtrans (const cs_di *A, int seed) ;
CXSPARSE_EXPORT int *cs_di_post (const int *parent, int n) ;
CXSPARSE_EXPORT cs_did *cs_di_scc (cs_di *A) ;
CXSPARSE_EXPORT int cs_di_scatter (const cs_di *A, int j, double beta, int *w, double *x,
    int mark, cs_di *C, int nz) ;
CXSPARSE_EXPORT int cs_di_tdfs (int j, int k, int *head, const int *next, int *post,
    int *stack) ;
CXSPARSE_EXPORT int cs_di_leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
    int *ancestor, int *jleaf) ;
CXSPARSE_EXPORT int cs_di_reach (cs_di *G, const cs_di *B, int k, int *xi, const int *pinv) ;
CXSPARSE_EXPORT int cs_di_spsolve (cs_di *L, const cs_di *B, int k, int *xi, double *x,
    const int *pinv, int lo) ;
CXSPARSE_EXPORT int cs_di_ereach (const cs_di *A, int k, const int *parent, int *s, int *w) ;
CXSPARSE_EXPORT int *cs_di_randperm (int n, int seed) ;

/* utilities */
CXSPARSE_EXPORT cs_did *cs_di_dalloc (int m, int n) ;
CXSPARSE_EXPORT cs_di *cs_di_done (cs_di *C, void *w, void *x, int ok) ;
CXSPARSE_EXPORT int *cs_di_idone (int *p, cs_di *C, void *w, int ok) ;
CXSPARSE_EXPORT cs_din *cs_di_ndone (cs_din *N, cs_di *C, void *w, void *x, int ok) ;
CXSPARSE_EXPORT cs_did *cs_di_ddone (cs_did *D, cs_di *C, void *w, int ok) ;


/* -------------------------------------------------------------------------- */
/* double/cs_long_t version of CXSparse */
/* -------------------------------------------------------------------------- */

/* --- primary CSparse routines and data structures ------------------------- */

typedef struct cs_dl_sparse  /* matrix in compressed-column or triplet form */
{
    cs_long_t nzmax ; /* maximum number of entries */
    cs_long_t m ;     /* number of rows */
    cs_long_t n ;     /* number of columns */
    cs_long_t *p ;    /* column pointers (size n+1) or col indlces (size nzmax) */
    cs_long_t *i ;    /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    cs_long_t nz ;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs_dl ;

CXSPARSE_EXPORT cs_dl *cs_dl_add (const cs_dl *A, const cs_dl *B, double alpha, double beta) ;
CXSPARSE_EXPORT cs_long_t cs_dl_cholsol (cs_long_t order, const cs_dl *A, double *b) ;
CXSPARSE_EXPORT cs_long_t cs_dl_dupl (cs_dl *A) ;
CXSPARSE_EXPORT cs_long_t cs_dl_entry (cs_dl *T, cs_long_t i, cs_long_t j, double x) ;
CXSPARSE_EXPORT cs_long_t cs_dl_lusol (cs_long_t order, const cs_dl *A, double *b, double tol) ;
CXSPARSE_EXPORT cs_long_t cs_dl_gaxpy (const cs_dl *A, const double *x, double *y) ;
CXSPARSE_EXPORT cs_dl *cs_dl_multiply (const cs_dl *A, const cs_dl *B) ;
CXSPARSE_EXPORT cs_long_t cs_dl_qrsol (cs_long_t order, const cs_dl *A, double *b) ;
CXSPARSE_EXPORT cs_dl *cs_dl_transpose (const cs_dl *A, cs_long_t values) ;
CXSPARSE_EXPORT cs_dl *cs_dl_compress (const cs_dl *T) ;
CXSPARSE_EXPORT double cs_dl_norm (const cs_dl *A) ;
CXSPARSE_EXPORT cs_long_t cs_dl_print (const cs_dl *A, cs_long_t brief) ;
CXSPARSE_EXPORT cs_dl *cs_dl_load (FILE *f) ;

/* utilities */
CXSPARSE_EXPORT void *cs_dl_calloc (cs_long_t n, size_t size) ;
CXSPARSE_EXPORT void *cs_dl_free (void *p) ;
CXSPARSE_EXPORT void *cs_dl_realloc (void *p, cs_long_t n, size_t size, cs_long_t *ok) ;
CXSPARSE_EXPORT cs_dl *cs_dl_spalloc (cs_long_t m, cs_long_t n, cs_long_t nzmax, cs_long_t values,
    cs_long_t t) ;
CXSPARSE_EXPORT cs_dl *cs_dl_spfree (cs_dl *A) ;
CXSPARSE_EXPORT cs_long_t cs_dl_sprealloc (cs_dl *A, cs_long_t nzmax) ;
CXSPARSE_EXPORT void *cs_dl_malloc (cs_long_t n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

typedef struct cs_dl_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    cs_long_t *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    cs_long_t *q ;        /* fill-reducing column permutation for LU and QR */
    cs_long_t *parent ;   /* elimination tree for Cholesky and QR */
    cs_long_t *cp ;       /* column pointers for Cholesky, row counts for QR */
    cs_long_t *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    cs_long_t m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;        /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;        /* # entries in U for LU; in R for QR */
} cs_dls ;

typedef struct cs_dl_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs_dl *L ;      /* L for LU and Cholesky, V for QR */
    cs_dl *U ;      /* U for LU, r for QR, not used for Cholesky */
    cs_long_t *pinv ; /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} cs_dln ;

typedef struct cs_dl_dmperm_results    /* cs_dl_dmperm or cs_dl_scc output */
{
    cs_long_t *p ;    /* size m, row permutation */
    cs_long_t *q ;    /* size n, column permutation */
    cs_long_t *r ;    /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    cs_long_t *s ;    /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    cs_long_t nb ;    /* # of blocks in fine dmperm decomposition */
    cs_long_t rr [5] ;    /* coarse row decomposition */
    cs_long_t cc [5] ;    /* coarse column decomposition */
} cs_dld ;

CXSPARSE_EXPORT cs_long_t *cs_dl_amd (cs_long_t order, const cs_dl *A) ;
CXSPARSE_EXPORT cs_dln *cs_dl_chol (const cs_dl *A, const cs_dls *S) ;
CXSPARSE_EXPORT cs_dld *cs_dl_dmperm (const cs_dl *A, cs_long_t seed) ;
CXSPARSE_EXPORT cs_long_t cs_dl_droptol (cs_dl *A, double tol) ;
CXSPARSE_EXPORT cs_long_t cs_dl_dropzeros (cs_dl *A) ;
CXSPARSE_EXPORT cs_long_t cs_dl_happly (const cs_dl *V, cs_long_t i, double beta, double *x) ;
CXSPARSE_EXPORT cs_long_t cs_dl_ipvec (const cs_long_t *p, const double *b, double *x, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t cs_dl_lsolve (const cs_dl *L, double *x) ;
CXSPARSE_EXPORT cs_long_t cs_dl_ltsolve (const cs_dl *L, double *x) ;
CXSPARSE_EXPORT cs_dln *cs_dl_lu (const cs_dl *A, const cs_dls *S, double tol) ;
CXSPARSE_EXPORT cs_dl *cs_dl_permute (const cs_dl *A, const cs_long_t *pinv, const cs_long_t *q,
    cs_long_t values) ;
CXSPARSE_EXPORT cs_long_t *cs_dl_pinv (const cs_long_t *p, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t cs_dl_pvec (const cs_long_t *p, const double *b, double *x, cs_long_t n) ;
CXSPARSE_EXPORT cs_dln *cs_dl_qr (const cs_dl *A, const cs_dls *S) ;
CXSPARSE_EXPORT cs_dls *cs_dl_schol (cs_long_t order, const cs_dl *A) ;
CXSPARSE_EXPORT cs_dls *cs_dl_sqr (cs_long_t order, const cs_dl *A, cs_long_t qr) ;
CXSPARSE_EXPORT cs_dl *cs_dl_symperm (const cs_dl *A, const cs_long_t *pinv, cs_long_t values) ;
CXSPARSE_EXPORT cs_long_t cs_dl_usolve (const cs_dl *U, double *x) ;
CXSPARSE_EXPORT cs_long_t cs_dl_utsolve (const cs_dl *U, double *x) ;
CXSPARSE_EXPORT cs_long_t cs_dl_updown (cs_dl *L, cs_long_t sigma, const cs_dl *C,
    const cs_long_t *parent) ;

/* utilities */
CXSPARSE_EXPORT cs_dls *cs_dl_sfree (cs_dls *S) ;
CXSPARSE_EXPORT cs_dln *cs_dl_nfree (cs_dln *N) ;
CXSPARSE_EXPORT cs_dld *cs_dl_dfree (cs_dld *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */

CXSPARSE_EXPORT cs_long_t *cs_dl_counts (const cs_dl *A, const cs_long_t *parent,
    const cs_long_t *post, cs_long_t ata) ;
CXSPARSE_EXPORT double cs_dl_cumsum (cs_long_t *p, cs_long_t *c, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t cs_dl_dfs (cs_long_t j, cs_dl *G, cs_long_t top, cs_long_t *xi,
    cs_long_t *pstack, const cs_long_t *pinv) ;
CXSPARSE_EXPORT cs_long_t *cs_dl_etree (const cs_dl *A, cs_long_t ata) ;
CXSPARSE_EXPORT cs_long_t cs_dl_fkeep (cs_dl *A,
    cs_long_t (*fkeep) (cs_long_t, cs_long_t, double, void *), void *other) ;
CXSPARSE_EXPORT double cs_dl_house (double *x, double *beta, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t *cs_dl_maxtrans (const cs_dl *A, cs_long_t seed) ;
CXSPARSE_EXPORT cs_long_t *cs_dl_post (const cs_long_t *parent, cs_long_t n) ;
CXSPARSE_EXPORT cs_dld *cs_dl_scc (cs_dl *A) ;
CXSPARSE_EXPORT cs_long_t cs_dl_scatter (const cs_dl *A, cs_long_t j, double beta, cs_long_t *w,
    double *x, cs_long_t mark,cs_dl *C, cs_long_t nz) ;
CXSPARSE_EXPORT cs_long_t cs_dl_tdfs (cs_long_t j, cs_long_t k, cs_long_t *head, const cs_long_t *next,
    cs_long_t *post, cs_long_t *stack) ;
CXSPARSE_EXPORT cs_long_t cs_dl_leaf (cs_long_t i, cs_long_t j, const cs_long_t *first,
    cs_long_t *maxfirst, cs_long_t *prevleaf, cs_long_t *ancestor, cs_long_t *jleaf) ;
CXSPARSE_EXPORT cs_long_t cs_dl_reach (cs_dl *G, const cs_dl *B, cs_long_t k, cs_long_t *xi,
    const cs_long_t *pinv) ;
CXSPARSE_EXPORT cs_long_t cs_dl_spsolve (cs_dl *L, const cs_dl *B, cs_long_t k, cs_long_t *xi,
    double *x, const cs_long_t *pinv, cs_long_t lo) ;
CXSPARSE_EXPORT cs_long_t cs_dl_ereach (const cs_dl *A, cs_long_t k, const cs_long_t *parent,
    cs_long_t *s, cs_long_t *w) ;
CXSPARSE_EXPORT cs_long_t *cs_dl_randperm (cs_long_t n, cs_long_t seed) ;

/* utilities */
CXSPARSE_EXPORT cs_dld *cs_dl_dalloc (cs_long_t m, cs_long_t n) ;
CXSPARSE_EXPORT cs_dl *cs_dl_done (cs_dl *C, void *w, void *x, cs_long_t ok) ;
CXSPARSE_EXPORT cs_long_t *cs_dl_idone (cs_long_t *p, cs_dl *C, void *w, cs_long_t ok) ;
CXSPARSE_EXPORT cs_dln *cs_dl_ndone (cs_dln *N, cs_dl *C, void *w, void *x, cs_long_t ok) ;
CXSPARSE_EXPORT cs_dld *cs_dl_ddone (cs_dld *D, cs_dl *C, void *w, cs_long_t ok) ;


/* -------------------------------------------------------------------------- */
/* complex/int version of CXSparse */
/* -------------------------------------------------------------------------- */

#ifndef NCOMPLEX

/* --- primary CSparse routines and data structures ------------------------- */

typedef struct cs_ci_sparse  /* matrix in compressed-column or triplet form */
{
    int nzmax ;     /* maximum number of entries */
    int m ;         /* number of rows */
    int n ;         /* number of columns */
    int *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;        /* row indices, size nzmax */
    cs_complex_t *x ;    /* numerical values, size nzmax */
    int nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs_ci ;

CXSPARSE_EXPORT cs_ci *cs_ci_add (const cs_ci *A, const cs_ci *B, cs_complex_t alpha,
    cs_complex_t beta) ;
CXSPARSE_EXPORT int cs_ci_cholsol (int order, const cs_ci *A, cs_complex_t *b) ;
CXSPARSE_EXPORT int cs_ci_dupl (cs_ci *A) ;
CXSPARSE_EXPORT int cs_ci_entry (cs_ci *T, int i, int j, cs_complex_t x) ;
CXSPARSE_EXPORT int cs_ci_lusol (int order, const cs_ci *A, cs_complex_t *b, double tol) ;
CXSPARSE_EXPORT int cs_ci_gaxpy (const cs_ci *A, const cs_complex_t *x, cs_complex_t *y) ;
CXSPARSE_EXPORT cs_ci *cs_ci_multiply (const cs_ci *A, const cs_ci *B) ;
CXSPARSE_EXPORT int cs_ci_qrsol (int order, const cs_ci *A, cs_complex_t *b) ;
CXSPARSE_EXPORT cs_ci *cs_ci_transpose (const cs_ci *A, int values) ;
CXSPARSE_EXPORT cs_ci *cs_ci_compress (const cs_ci *T) ;
CXSPARSE_EXPORT double cs_ci_norm (const cs_ci *A) ;
CXSPARSE_EXPORT int cs_ci_print (const cs_ci *A, int brief) ;
CXSPARSE_EXPORT cs_ci *cs_ci_load (FILE *f) ;

/* utilities */
CXSPARSE_EXPORT void *cs_ci_calloc (int n, size_t size) ;
CXSPARSE_EXPORT void *cs_ci_free (void *p) ;
CXSPARSE_EXPORT void *cs_ci_realloc (void *p, int n, size_t size, int *ok) ;
CXSPARSE_EXPORT cs_ci *cs_ci_spalloc (int m, int n, int nzmax, int values, int t) ;
CXSPARSE_EXPORT cs_ci *cs_ci_spfree (cs_ci *A) ;
CXSPARSE_EXPORT int cs_ci_sprealloc (cs_ci *A, int nzmax) ;
CXSPARSE_EXPORT void *cs_ci_malloc (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

typedef struct cs_ci_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    int *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    int *q ;        /* fill-reducing column permutation for LU and QR */
    int *parent ;   /* elimination tree for Cholesky and QR */
    int *cp ;       /* column pointers for Cholesky, row counts for QR */
    int *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    int m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;    /* # entries in U for LU; in R for QR */
} cs_cis ;

typedef struct cs_ci_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs_ci *L ;      /* L for LU and Cholesky, V for QR */
    cs_ci *U ;      /* U for LU, r for QR, not used for Cholesky */
    int *pinv ;     /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} cs_cin ;

typedef struct cs_ci_dmperm_results    /* cs_ci_dmperm or cs_ci_scc output */
{
    int *p ;        /* size m, row permutation */
    int *q ;        /* size n, column permutation */
    int *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    int *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    int nb ;        /* # of blocks in fine dmperm decomposition */
    int rr [5] ;    /* coarse row decomposition */
    int cc [5] ;    /* coarse column decomposition */
} cs_cid ;

CXSPARSE_EXPORT int *cs_ci_amd (int order, const cs_ci *A) ;
CXSPARSE_EXPORT cs_cin *cs_ci_chol (const cs_ci *A, const cs_cis *S) ;
CXSPARSE_EXPORT cs_cid *cs_ci_dmperm (const cs_ci *A, int seed) ;
CXSPARSE_EXPORT int cs_ci_droptol (cs_ci *A, double tol) ;
CXSPARSE_EXPORT int cs_ci_dropzeros (cs_ci *A) ;
CXSPARSE_EXPORT int cs_ci_happly (const cs_ci *V, int i, double beta, cs_complex_t *x) ;
CXSPARSE_EXPORT int cs_ci_ipvec (const int *p, const cs_complex_t *b, cs_complex_t *x, int n) ;
CXSPARSE_EXPORT int cs_ci_lsolve (const cs_ci *L, cs_complex_t *x) ;
CXSPARSE_EXPORT int cs_ci_ltsolve (const cs_ci *L, cs_complex_t *x) ;
CXSPARSE_EXPORT cs_cin *cs_ci_lu (const cs_ci *A, const cs_cis *S, double tol) ;
CXSPARSE_EXPORT cs_ci *cs_ci_permute (const cs_ci *A, const int *pinv, const int *q,
    int values) ;
CXSPARSE_EXPORT int *cs_ci_pinv (const int *p, int n) ;
CXSPARSE_EXPORT int cs_ci_pvec (const int *p, const cs_complex_t *b, cs_complex_t *x, int n) ;
CXSPARSE_EXPORT cs_cin *cs_ci_qr (const cs_ci *A, const cs_cis *S) ;
CXSPARSE_EXPORT cs_cis *cs_ci_schol (int order, const cs_ci *A) ;
CXSPARSE_EXPORT cs_cis *cs_ci_sqr (int order, const cs_ci *A, int qr) ;
CXSPARSE_EXPORT cs_ci *cs_ci_symperm (const cs_ci *A, const int *pinv, int values) ;
CXSPARSE_EXPORT int cs_ci_usolve (const cs_ci *U, cs_complex_t *x) ;
CXSPARSE_EXPORT int cs_ci_utsolve (const cs_ci *U, cs_complex_t *x) ;
CXSPARSE_EXPORT int cs_ci_updown (cs_ci *L, int sigma, const cs_ci *C, const int *parent) ;

/* utilities */
CXSPARSE_EXPORT cs_cis *cs_ci_sfree (cs_cis *S) ;
CXSPARSE_EXPORT cs_cin *cs_ci_nfree (cs_cin *N) ;
CXSPARSE_EXPORT cs_cid *cs_ci_dfree (cs_cid *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */

CXSPARSE_EXPORT int *cs_ci_counts (const cs_ci *A, const int *parent, const int *post,
    int ata) ;
CXSPARSE_EXPORT double cs_ci_cumsum (int *p, int *c, int n) ;
CXSPARSE_EXPORT int cs_ci_dfs (int j, cs_ci *G, int top, int *xi, int *pstack,
    const int *pinv) ;
CXSPARSE_EXPORT int *cs_ci_etree (const cs_ci *A, int ata) ;
CXSPARSE_EXPORT int cs_ci_fkeep (cs_ci *A, int (*fkeep) (int, int, cs_complex_t, void *),
    void *other) ;
CXSPARSE_EXPORT cs_complex_t cs_ci_house (cs_complex_t *x, double *beta, int n) ;
CXSPARSE_EXPORT int *cs_ci_maxtrans (const cs_ci *A, int seed) ;
CXSPARSE_EXPORT int *cs_ci_post (const int *parent, int n) ;
CXSPARSE_EXPORT cs_cid *cs_ci_scc (cs_ci *A) ;
CXSPARSE_EXPORT int cs_ci_scatter (const cs_ci *A, int j, cs_complex_t beta, int *w,
    cs_complex_t *x, int mark,cs_ci *C, int nz) ;
CXSPARSE_EXPORT int cs_ci_tdfs (int j, int k, int *head, const int *next, int *post,
    int *stack) ;
CXSPARSE_EXPORT int cs_ci_leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
    int *ancestor, int *jleaf) ;
CXSPARSE_EXPORT int cs_ci_reach (cs_ci *G, const cs_ci *B, int k, int *xi, const int *pinv) ;
CXSPARSE_EXPORT int cs_ci_spsolve (cs_ci *L, const cs_ci *B, int k, int *xi,
    cs_complex_t *x, const int *pinv, int lo) ;
CXSPARSE_EXPORT int cs_ci_ereach (const cs_ci *A, int k, const int *parent, int *s, int *w) ;
CXSPARSE_EXPORT int *cs_ci_randperm (int n, int seed) ;

/* utilities */
CXSPARSE_EXPORT cs_cid *cs_ci_dalloc (int m, int n) ;
CXSPARSE_EXPORT cs_ci *cs_ci_done (cs_ci *C, void *w, void *x, int ok) ;
CXSPARSE_EXPORT int *cs_ci_idone (int *p, cs_ci *C, void *w, int ok) ;
CXSPARSE_EXPORT cs_cin *cs_ci_ndone (cs_cin *N, cs_ci *C, void *w, void *x, int ok) ;
CXSPARSE_EXPORT cs_cid *cs_ci_ddone (cs_cid *D, cs_ci *C, void *w, int ok) ;


/* -------------------------------------------------------------------------- */
/* complex/cs_long_t version of CXSparse */
/* -------------------------------------------------------------------------- */

/* --- primary CSparse routines and data structures ------------------------- */

typedef struct cs_cl_sparse  /* matrix in compressed-column or triplet form */
{
    cs_long_t nzmax ; /* maximum number of entries */
    cs_long_t m ;     /* number of rows */
    cs_long_t n ;     /* number of columns */
    cs_long_t *p ;    /* column pointers (size n+1) or col indlces (size nzmax) */
    cs_long_t *i ;    /* row indices, size nzmax */
    cs_complex_t *x ;    /* numerical values, size nzmax */
    cs_long_t nz ;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs_cl ;

CXSPARSE_EXPORT cs_cl *cs_cl_add (const cs_cl *A, const cs_cl *B, cs_complex_t alpha,
    cs_complex_t beta) ;
CXSPARSE_EXPORT cs_long_t cs_cl_cholsol (cs_long_t order, const cs_cl *A, cs_complex_t *b) ;
CXSPARSE_EXPORT cs_long_t cs_cl_dupl (cs_cl *A) ;
CXSPARSE_EXPORT cs_long_t cs_cl_entry (cs_cl *T, cs_long_t i, cs_long_t j, cs_complex_t x) ;
CXSPARSE_EXPORT cs_long_t cs_cl_lusol (cs_long_t order, const cs_cl *A, cs_complex_t *b,
    double tol) ;
CXSPARSE_EXPORT cs_long_t cs_cl_gaxpy (const cs_cl *A, const cs_complex_t *x, cs_complex_t *y) ;
CXSPARSE_EXPORT cs_cl *cs_cl_multiply (const cs_cl *A, const cs_cl *B) ;
CXSPARSE_EXPORT cs_long_t cs_cl_qrsol (cs_long_t order, const cs_cl *A, cs_complex_t *b) ;
CXSPARSE_EXPORT cs_cl *cs_cl_transpose (const cs_cl *A, cs_long_t values) ;
CXSPARSE_EXPORT cs_cl *cs_cl_compress (const cs_cl *T) ;
CXSPARSE_EXPORT double cs_cl_norm (const cs_cl *A) ;
CXSPARSE_EXPORT cs_long_t cs_cl_print (const cs_cl *A, cs_long_t brief) ;
CXSPARSE_EXPORT cs_cl *cs_cl_load (FILE *f) ;

/* utilities */
CXSPARSE_EXPORT void *cs_cl_calloc (cs_long_t n, size_t size) ;
CXSPARSE_EXPORT void *cs_cl_free (void *p) ;
CXSPARSE_EXPORT void *cs_cl_realloc (void *p, cs_long_t n, size_t size, cs_long_t *ok) ;
CXSPARSE_EXPORT cs_cl *cs_cl_spalloc (cs_long_t m, cs_long_t n, cs_long_t nzmax, cs_long_t values,
    cs_long_t t) ;
CXSPARSE_EXPORT cs_cl *cs_cl_spfree (cs_cl *A) ;
CXSPARSE_EXPORT cs_long_t cs_cl_sprealloc (cs_cl *A, cs_long_t nzmax) ;
CXSPARSE_EXPORT void *cs_cl_malloc (cs_long_t n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

typedef struct cs_cl_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    cs_long_t *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    cs_long_t *q ;        /* fill-reducing column permutation for LU and QR */
    cs_long_t *parent ;   /* elimination tree for Cholesky and QR */
    cs_long_t *cp ;       /* column pointers for Cholesky, row counts for QR */
    cs_long_t *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    cs_long_t m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;        /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;        /* # entries in U for LU; in R for QR */
} cs_cls ;

typedef struct cs_cl_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs_cl *L ;          /* L for LU and Cholesky, V for QR */
    cs_cl *U ;          /* U for LU, r for QR, not used for Cholesky */
    cs_long_t *pinv ;     /* partial pivoting for LU */
    double *B ;         /* beta [0..n-1] for QR */
} cs_cln ;

typedef struct cs_cl_dmperm_results    /* cs_cl_dmperm or cs_cl_scc output */
{
    cs_long_t *p ;    /* size m, row permutation */
    cs_long_t *q ;    /* size n, column permutation */
    cs_long_t *r ;    /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    cs_long_t *s ;    /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    cs_long_t nb ;    /* # of blocks in fine dmperm decomposition */
    cs_long_t rr [5] ;   /* coarse row decomposition */
    cs_long_t cc [5] ;   /* coarse column decomposition */
} cs_cld ;

CXSPARSE_EXPORT cs_long_t *cs_cl_amd (cs_long_t order, const cs_cl *A) ;
CXSPARSE_EXPORT cs_cln *cs_cl_chol (const cs_cl *A, const cs_cls *S) ;
CXSPARSE_EXPORT cs_cld *cs_cl_dmperm (const cs_cl *A, cs_long_t seed) ;
CXSPARSE_EXPORT cs_long_t cs_cl_droptol (cs_cl *A, double tol) ;
CXSPARSE_EXPORT cs_long_t cs_cl_dropzeros (cs_cl *A) ;
CXSPARSE_EXPORT cs_long_t cs_cl_happly (const cs_cl *V, cs_long_t i, double beta, cs_complex_t *x) ;
CXSPARSE_EXPORT cs_long_t cs_cl_ipvec (const cs_long_t *p, const cs_complex_t *b,
    cs_complex_t *x, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t cs_cl_lsolve (const cs_cl *L, cs_complex_t *x) ;
CXSPARSE_EXPORT cs_long_t cs_cl_ltsolve (const cs_cl *L, cs_complex_t *x) ;
CXSPARSE_EXPORT cs_cln *cs_cl_lu (const cs_cl *A, const cs_cls *S, double tol) ;
CXSPARSE_EXPORT cs_cl *cs_cl_permute (const cs_cl *A, const cs_long_t *pinv, const cs_long_t *q,
    cs_long_t values) ;
CXSPARSE_EXPORT cs_long_t *cs_cl_pinv (const cs_long_t *p, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t cs_cl_pvec (const cs_long_t *p, const cs_complex_t *b,
    cs_complex_t *x, cs_long_t n) ;
CXSPARSE_EXPORT cs_cln *cs_cl_qr (const cs_cl *A, const cs_cls *S) ;
CXSPARSE_EXPORT cs_cls *cs_cl_schol (cs_long_t order, const cs_cl *A) ;
CXSPARSE_EXPORT cs_cls *cs_cl_sqr (cs_long_t order, const cs_cl *A, cs_long_t qr) ;
CXSPARSE_EXPORT cs_cl *cs_cl_symperm (const cs_cl *A, const cs_long_t *pinv, cs_long_t values) ;
CXSPARSE_EXPORT cs_long_t cs_cl_usolve (const cs_cl *U, cs_complex_t *x) ;
CXSPARSE_EXPORT cs_long_t cs_cl_utsolve (const cs_cl *U, cs_complex_t *x) ;
CXSPARSE_EXPORT cs_long_t cs_cl_updown (cs_cl *L, cs_long_t sigma, const cs_cl *C,
    const cs_long_t *parent) ;

/* utilities */
CXSPARSE_EXPORT cs_cls *cs_cl_sfree (cs_cls *S) ;
CXSPARSE_EXPORT cs_cln *cs_cl_nfree (cs_cln *N) ;
CXSPARSE_EXPORT cs_cld *cs_cl_dfree (cs_cld *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */

CXSPARSE_EXPORT cs_long_t *cs_cl_counts (const cs_cl *A, const cs_long_t *parent,
    const cs_long_t *post, cs_long_t ata) ;
CXSPARSE_EXPORT double cs_cl_cumsum (cs_long_t *p, cs_long_t *c, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t cs_cl_dfs (cs_long_t j, cs_cl *G, cs_long_t top, cs_long_t *xi,
    cs_long_t *pstack, const cs_long_t *pinv) ;
CXSPARSE_EXPORT cs_long_t *cs_cl_etree (const cs_cl *A, cs_long_t ata) ;
CXSPARSE_EXPORT cs_long_t cs_cl_fkeep (cs_cl *A,
    cs_long_t (*fkeep) (cs_long_t, cs_long_t, cs_complex_t, void *), void *other) ;
CXSPARSE_EXPORT cs_complex_t cs_cl_house (cs_complex_t *x, double *beta, cs_long_t n) ;
CXSPARSE_EXPORT cs_long_t *cs_cl_maxtrans (const cs_cl *A, cs_long_t seed) ;
CXSPARSE_EXPORT cs_long_t *cs_cl_post (const cs_long_t *parent, cs_long_t n) ;
CXSPARSE_EXPORT cs_cld *cs_cl_scc (cs_cl *A) ;
CXSPARSE_EXPORT cs_long_t cs_cl_scatter (const cs_cl *A, cs_long_t j, cs_complex_t beta,
    cs_long_t *w, cs_complex_t *x, cs_long_t mark,cs_cl *C, cs_long_t nz) ;
CXSPARSE_EXPORT cs_long_t cs_cl_tdfs (cs_long_t j, cs_long_t k, cs_long_t *head, const cs_long_t *next,
    cs_long_t *post, cs_long_t *stack) ;
CXSPARSE_EXPORT cs_long_t cs_cl_leaf (cs_long_t i, cs_long_t j, const cs_long_t *first,
    cs_long_t *maxfirst, cs_long_t *prevleaf, cs_long_t *ancestor, cs_long_t *jleaf) ;
CXSPARSE_EXPORT cs_long_t cs_cl_reach (cs_cl *G, const cs_cl *B, cs_long_t k, cs_long_t *xi,
    const cs_long_t *pinv) ;
CXSPARSE_EXPORT cs_long_t cs_cl_spsolve (cs_cl *L, const cs_cl *B, cs_long_t k, cs_long_t *xi,
    cs_complex_t *x, const cs_long_t *pinv, cs_long_t lo) ;
CXSPARSE_EXPORT cs_long_t cs_cl_ereach (const cs_cl *A, cs_long_t k, const cs_long_t *parent,
    cs_long_t *s, cs_long_t *w) ;
CXSPARSE_EXPORT cs_long_t *cs_cl_randperm (cs_long_t n, cs_long_t seed) ;

/* utilities */
CXSPARSE_EXPORT cs_cld *cs_cl_dalloc (cs_long_t m, cs_long_t n) ;
CXSPARSE_EXPORT cs_cl *cs_cl_done (cs_cl *C, void *w, void *x, cs_long_t ok) ;
CXSPARSE_EXPORT cs_long_t *cs_cl_idone (cs_long_t *p, cs_cl *C, void *w, cs_long_t ok) ;
CXSPARSE_EXPORT cs_cln *cs_cl_ndone (cs_cln *N, cs_cl *C, void *w, void *x, cs_long_t ok) ;
CXSPARSE_EXPORT cs_cld *cs_cl_ddone (cs_cld *D, cs_cl *C, void *w, cs_long_t ok) ;

#endif

/* -------------------------------------------------------------------------- */
/* Macros for constructing each version of CSparse */
/* -------------------------------------------------------------------------- */

#ifdef CS_LONG
#define CS_INT cs_long_t
#define CS_INT_MAX cs_long_t_max
#define CS_ID cs_long_t_id
#ifdef CS_COMPLEX
#define CS_ENTRY cs_complex_t
#define CS_NAME(nm) cs_cl ## nm
#define cs cs_cl
#else
#define CS_ENTRY double
#define CS_NAME(nm) cs_dl ## nm
#define cs cs_dl
#endif
#else
#define CS_INT int
#define CS_INT_MAX INT_MAX
#define CS_ID "%d"
#ifdef CS_COMPLEX
#define CS_ENTRY cs_complex_t
#define CS_NAME(nm) cs_ci ## nm
#define cs cs_ci
#else
#define CS_ENTRY double
#define CS_NAME(nm) cs_di ## nm
#define cs cs_di
#endif
#endif

#ifdef CS_COMPLEX
#ifndef NCOMPLEX
#define CS_REAL(x) creal(x)
#define CS_IMAG(x) cimag(x)
#define CS_CONJ(x) conj(x)
#define CS_ABS(x) cabs(x)
#elif defined(__cplusplus)
#define CS_REAL(x) (x).real()
#define CS_IMAG(x) (x).imag()
#define CS_CONJ(x) conj(x)
#define CS_ABS(x) abs(x)
#endif
#else
#define CS_REAL(x) (x)
#define CS_IMAG(x) (0.)
#define CS_CONJ(x) (x)
#define CS_ABS(x) fabs(x)
#endif

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

/* --- primary CSparse routines and data structures ------------------------- */

#define cs_add CS_NAME (_add)
#define cs_cholsol CS_NAME (_cholsol)
#define cs_dupl CS_NAME (_dupl)
#define cs_entry CS_NAME (_entry)
#define cs_lusol CS_NAME (_lusol)
#define cs_gaxpy CS_NAME (_gaxpy)
#define cs_multiply CS_NAME (_multiply)
#define cs_qrsol CS_NAME (_qrsol)
#define cs_transpose CS_NAME (_transpose)
#define cs_compress CS_NAME (_compress)
#define cs_norm CS_NAME (_norm)
#define cs_print CS_NAME (_print)
#define cs_load CS_NAME (_load)

/* utilities */
#define cs_calloc CS_NAME (_calloc)
#define cs_free CS_NAME (_free)
#define cs_realloc CS_NAME (_realloc)
#define cs_spalloc CS_NAME (_spalloc)
#define cs_spfree CS_NAME (_spfree)
#define cs_sprealloc CS_NAME (_sprealloc)
#define cs_malloc CS_NAME (_malloc)

/* --- secondary CSparse routines and data structures ----------------------- */
#define css CS_NAME (s)
#define csn CS_NAME (n)
#define csd CS_NAME (d)

#define cs_amd CS_NAME (_amd)
#define cs_chol CS_NAME (_chol)
#define cs_dmperm CS_NAME (_dmperm)
#define cs_droptol CS_NAME (_droptol)
#define cs_dropzeros CS_NAME (_dropzeros)
#define cs_happly CS_NAME (_happly)
#define cs_ipvec CS_NAME (_ipvec)
#define cs_lsolve CS_NAME (_lsolve)
#define cs_ltsolve CS_NAME (_ltsolve)
#define cs_lu CS_NAME (_lu)
#define cs_permute CS_NAME (_permute)
#define cs_pinv CS_NAME (_pinv)
#define cs_pvec CS_NAME (_pvec)
#define cs_qr CS_NAME (_qr)
#define cs_schol CS_NAME (_schol)
#define cs_sqr CS_NAME (_sqr)
#define cs_symperm CS_NAME (_symperm)
#define cs_usolve CS_NAME (_usolve)
#define cs_utsolve CS_NAME (_utsolve)
#define cs_updown CS_NAME (_updown)

/* utilities */
#define cs_sfree CS_NAME (_sfree)
#define cs_nfree CS_NAME (_nfree)
#define cs_dfree CS_NAME (_dfree)

/* --- tertiary CSparse routines -------------------------------------------- */
#define cs_counts CS_NAME (_counts)
#define cs_cumsum CS_NAME (_cumsum)
#define cs_dfs CS_NAME (_dfs)
#define cs_etree CS_NAME (_etree)
#define cs_fkeep CS_NAME (_fkeep)
#define cs_house CS_NAME (_house)
#define cs_invmatch CS_NAME (_invmatch)
#define cs_maxtrans CS_NAME (_maxtrans)
#define cs_post CS_NAME (_post)
#define cs_scc CS_NAME (_scc)
#define cs_scatter CS_NAME (_scatter)
#define cs_tdfs CS_NAME (_tdfs)
#define cs_reach CS_NAME (_reach)
#define cs_spsolve CS_NAME (_spsolve)
#define cs_ereach CS_NAME (_ereach)
#define cs_randperm CS_NAME (_randperm)
#define cs_leaf CS_NAME (_leaf)

/* utilities */
#define cs_dalloc CS_NAME (_dalloc)
#define cs_done CS_NAME (_done)
#define cs_idone CS_NAME (_idone)
#define cs_ndone CS_NAME (_ndone)
#define cs_ddone CS_NAME (_ddone)

/* -------------------------------------------------------------------------- */
/* Conversion routines */
/* -------------------------------------------------------------------------- */

#ifndef NCOMPLEX
CXSPARSE_EXPORT cs_di *cs_i_real (cs_ci *A, int real) ;
CXSPARSE_EXPORT cs_ci *cs_i_complex (cs_di *A, int real) ;
CXSPARSE_EXPORT cs_dl *cs_l_real (cs_cl *A, cs_long_t real) ;
CXSPARSE_EXPORT cs_cl *cs_l_complex (cs_dl *A, cs_long_t real) ;
#endif

#if defined(__cplusplus)
}
#endif
#endif
