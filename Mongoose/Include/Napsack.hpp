#ifndef NAPSACK_HPP_
#define NAPSACK_HPP_

#include "Mongoose_internal.hpp"

namespace SuiteSparse_Mongoose
{

void QPmaxheap_build
(
    Int *heap,                      /* on input, an unsorted set of elements */
    Int size,                       /* size of the heap */
    Double *x
);

Int QPmaxheap_delete             /* return new size of heap */
(
    Int *heap,                   /* containing indices into x, 1..n on input */
    Int size,                    /* size of the heap */
    Double *x                    /* not modified */
);

void QPmaxheapify
(
    Int p,                       /* start at node p in the heap */
    Int *heap,                   /* size n, containing indices into x */
    Int size,                    /* heap [ ... nheap] is in use */
    Double *x                    /* not modified */
);

Int QPmaxheap_add
(
    Int leaf ,   /* the new leaf */
    Int *heap ,  /* size n, containing indices into x */
    Double *x ,  /* not modified */
    Int size     /* number of elements in heap not counting new one */
);

void QPmaxheap_check
(
    Int *heap,  /* vector of size n+1 */
    Double *x,  /* vector of size n */
    Int size,       /* # items in heap */
    Int n,
    Int p       /* start checking at heap [p] */
);

void QPminheap_build
(
    Int *heap,                      /* on input, an unsorted set of elements */
    Int size,                       /* size of the heap */
    Double *x
);

Int QPminheap_delete             /* return new size of heap */
(
    Int *heap,                   /* containing indices into x, 1..n on input */
    Int size,                    /* size of the heap */
    Double *x                    /* not modified */
);

void QPminheapify
(
    Int p,                       /* start at node p in the heap */
    Int *heap,                   /* size n, containing indices into x */
    Int size,                    /* heap [ ... nheap] is in use */
    Double *x                    /* not modified */
);

Int QPminheap_add
(
    Int leaf ,   /* the new leaf */
    Int *heap ,  /* size n, containing indices into x */
    Double *x ,  /* not modified */
    Int size     /* number of elements in heap not counting new one */
);

void QPminheap_check
(
    Int *heap,  /* vector of size n+1 */
    Double *x,  /* vector of size n */
    Int size,       /* # items in heap */
    Int n,
    Int p       /* start checking at heap [p] */
);

Double QPnapsack            /* return the final lambda */
(
    Double *x,              /* holds y on input, and the solution x on output */
    Int n,                  /* size of x, constraint lo <= a'x <= hi */
    Double lo,              /* partition lower bound */
    Double hi,              /* partition upper bound */
    Double *a,              /* vector of nodal weights */
    Double Lambda,          /* initial guess for lambda */
    Int *ix,                /* ix_i = +1,-1, or 0 on input, x_i =1,0, or 0< x_i< 1*/
    Double *w,              /* work array of size n */
    Int *heap1,             /* work array of size n+1 */
    Int *heap2              /* work array of size n+1 */
);

void QPnapsack2
(
    Double *x,         /* holds y on input, and the solution x on output */
    Double *lambda,    /* initial guess (input) final value (output) multiplier
                          lambda [0], lambda [1] */
    Int *ix,           /* ix_i = +1, -1, or 0 on input, x_i =1, 0, or 0 < x_i < 1 */
    Int n,             /* # cols ? */
    Int m,             /* # rows ? */
    Double *a,         /* edge weights ? */
    Int *lo,           /* Problem lo limits in 0, 1 */
    Int *hi,           /* Problem hi limits in 0, 1 */
    Double *wx1,       /* Double workspace #1 */
    Int *wi1,          /* Integer workspace #1 */
    Int *wi2           /* Integer workspace #2 */
);

Double QPnapup              /* return lambda */
(
    Double *x,              /* holds y on input, not modified */
    Int n,                  /* size of x */
    Double lambda,          /* initial guess for the shift */
    Double *a,              /* input constraint vector */
    Double b,               /* input constraint scalar */
    Double *breakpts,       /* break points */
    Int *bound_heap,        /* work array */
    Int *free_heap          /* work array */
);

Double QPnapdown            /* return lambda */
(
    Double *x,              /* holds y on input, not modified */
    Int n,                  /* size of x */
    Double lambda,          /* initial guess for the shift */
    Double *a,              /* input constraint vector */
    Double b,               /* input constraint scalar */
    Double *breakpts,       /* break points */
    Int *bound_heap,        /* work array */
    Int *free_heap          /* work array */
);

}

#endif
