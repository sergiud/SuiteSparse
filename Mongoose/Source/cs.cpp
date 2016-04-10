
#include "Mongoose_internal.hpp"
#include "cs.hpp"
#include "SuiteSparse_config.h"

namespace SuiteSparse_Mongoose
{

csi *cs_counts (const cs *A, const csi *parent, const csi *post, csi ata) ;
double cs_cumsum (csi *p, csi *c, csi n) ;
csi cs_scatter (const cs *A, csi j, double beta, csi *w, double *x, csi mark, cs *C, csi nz);
csi cs_sprealloc (cs *A, csi nzmax) ;
void *cs_realloc (void *p, csi n, size_t size, csi *ok) ;
cs *cs_done (cs *C, void *w, void *x, csi ok) ;

//-----------------------------------------------------------------------------
// add an entry to a triplet matrix; return 1 if ok, 0 otherwise
//-----------------------------------------------------------------------------
csi cs_entry (cs *T, csi i, csi j, double x)
{
    if (!CS_TRIPLET (T) || i < 0 || j < 0) return (0) ;     /* check inputs */
    if (T->nz >= T->nzmax && !cs_sprealloc (T,2*(T->nzmax))) return (0) ;
    if (T->x) T->x [T->nz] = x ;
    T->i [T->nz] = i ;
    T->p [T->nz++] = j ;
    T->m = MONGOOSE_MAX2 (T->m, i+1) ;
    T->n = MONGOOSE_MAX2 (T->n, j+1) ;
    return (1) ;
}

//-----------------------------------------------------------------------------
// C = A'
//-----------------------------------------------------------------------------
cs *cs_transpose (
    const cs *A,
    csi values
)
{
    csi p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
    double *Cx, *Ax ;
    cs *C ;
    if (!CS_CSC (A)) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = (cs*) cs_spalloc (n, m, Ap [n], values && Ax, 0) ; /* allocate result */
    w = (csi*) SuiteSparse_calloc (m, sizeof (csi)) ;      /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;       /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    cs_cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return (cs_done (C, w, NULL, 1)) ;  /* success; release w and return C */
}

//-----------------------------------------------------------------------------
// C = alpha*A + beta*B
//-----------------------------------------------------------------------------
cs *cs_add (const cs *A, const cs *B, double alpha, double beta)
{
    csi p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
    double *x, *Bx, *Cx ;
    cs *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;         /* check inputs */
    if (A->m != B->m || A->n != B->n) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
    w = (csi*) SuiteSparse_calloc (m, sizeof (csi)) ;              /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? (double*) SuiteSparse_malloc (m, sizeof (double)) : NULL ;    /* get workspace */
    C = cs_spalloc (m, n, anz + bnz, values, 0) ;           /* allocate result*/
    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)
    {
        Cp [j] = nz ;                   /* column j of C starts here */
        nz = cs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        nz = cs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    return ((cs*) cs_done (C, w, x, 1)) ;     /* success; release workspace, return C */
}

//-----------------------------------------------------------------------------
// C = compressed-column form of a triplet matrix T
//-----------------------------------------------------------------------------
cs *cs_compress (const cs *T)
{
    csi m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    double *Cx, *Tx ;
    cs *C ;
    if (!CS_TRIPLET (T)) return (NULL) ;                /* check inputs */
    m = T->m ; n = T->n ; Ti = T->i ; Tj = T->p ; Tx = T->x ; nz = T->nz ;
    C = cs_spalloc (m, n, nz, Tx != NULL, 0) ;          /* allocate result */
    w = (csi*) SuiteSparse_calloc (n, sizeof (csi)) ;                   /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;           /* column counts */
    cs_cumsum (Cp, w, n) ;                              /* column pointers */
    for (k = 0 ; k < nz ; k++)
    {
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
    }
    return (cs_done (C, w, NULL, 1)) ;      /* success; release w and return C */
}

//-----------------------------------------------------------------------------
// p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c
//-----------------------------------------------------------------------------
double cs_cumsum (csi *p, csi *c, csi n)
{
    csi i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid csi overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

#if 0
//-----------------------------------------------------------------------------
// load a triplet matrix from a file
//-----------------------------------------------------------------------------
cs *cs_load (FILE *f)
{
    double i, j ;   /* use double for integers to avoid csi conflicts */
    double x ;
    cs *T ;
    if (!f) return (NULL) ;                             /* check inputs */
    T = cs_spalloc (0, 0, 1, 1, 1) ;                    /* allocate result */
    while (fscanf (f, "%lg %lg %lg\n", &i, &j, &x) == 3)
    {
        if (!cs_entry (T, (csi) i, (csi) j, x)) return (cs_spfree (T)) ;
    }
    return (T) ;
}
#endif

//-----------------------------------------------------------------------------
// x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse
//-----------------------------------------------------------------------------
csi cs_scatter (const cs *A, csi j, double beta, csi *w, double *x, csi mark,
    cs *C, csi nz)
{
    csi i, p, *Ap, *Ai, *Ci ;
    double *Ax ;
    if (!CS_CSC (A) || !w || !CS_CSC (C)) return (-1) ;     /* check inputs */
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}

//-----------------------------------------------------------------------------
// utility functions
//-----------------------------------------------------------------------------
/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet)
{
    cs *A = (cs*) SuiteSparse_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = MONGOOSE_MAX2 (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = (csi*) SuiteSparse_malloc (triplet ? nzmax : n+1, sizeof (csi)) ;
    A->i = (csi*) SuiteSparse_malloc (nzmax, sizeof (csi)) ;
    A->x = values ? (double*) SuiteSparse_malloc (nzmax, sizeof (double)) : NULL ;
    return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
csi cs_sprealloc (cs *A, csi nzmax)
{
    csi ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (CS_CSC (A)) ? (A->p [A->n]) : A->nz ;
    A->i = (csi*) cs_realloc (A->i, nzmax, sizeof (csi), &oki) ;
    if (CS_TRIPLET (A)) A->p = (csi*) cs_realloc (A->p, nzmax, sizeof (csi), &okj) ;
    if (A->x) A->x = (double*) cs_realloc (A->x, nzmax, sizeof (double), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* release a sparse matrix */
cs *cs_spfree (cs *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    SuiteSparse_free (A->p) ;
    SuiteSparse_free (A->i) ;
    SuiteSparse_free (A->x) ;
    return ((cs *) SuiteSparse_free (A)) ;   /* release the cs struct and return NULL */
}

/* wrapper for realloc */
void *cs_realloc (void *p, csi n, size_t size, csi *ok)
{
    void *pnew ;
    pnew = realloc (p, MONGOOSE_MAX2 (n,1) * size) ;/* realloc the block */
    *ok = (pnew != NULL) ;                    /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;               /* return original p if failure */
}

/* release workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, csi ok)
{
    SuiteSparse_free (w) ;                       /* release workspace */
    SuiteSparse_free (x) ;
    return (ok ? C : cs_spfree (C)) ;   /* return result if OK, else release it */
}

}