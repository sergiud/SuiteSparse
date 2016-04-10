#ifndef _CS_H
#define _CS_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#ifdef MATLAB_MEX_FILE
#undef csi
#define csi mwSignedIndex
#endif
#ifndef csi
#define csi ptrdiff_t
#endif

/* CSparse Macros */
#ifndef CS_CSC
#define CS_CSC(A) (A && (A->nz == -1))
#endif
#ifndef CS_TRIPLET
#define CS_TRIPLET(A) (A && (A->nz >= 0))
#endif

namespace SuiteSparse_Mongoose
{

/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    csi nzmax ;     /* maximum number of entries */
    csi m ;         /* number of rows */
    csi n ;         /* number of columns */
    csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    csi *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    csi nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
cs *cs_transpose (const cs *A, csi values) ;

#if 0
cs *cs_load (FILE *f) ;
#endif

csi cs_entry (cs *T, csi i, csi j, double x) ;
cs *cs_compress (const cs *T) ;

cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet) ;
cs *cs_spfree (cs *A) ;

}

#endif