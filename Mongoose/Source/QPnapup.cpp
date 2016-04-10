/* ========================================================================== */
/* === QPnapup ============================================================== */
/* ========================================================================== */

/* Find x that minimizes ||x-y|| while satisfying the constraints
   0 <= x <= 1, a'x = b. If a = NULL, then it is assumed that a = 1.
   The algorithm is described in the napsack comments.
   It is assumed that the starting guess lambda for the dual multiplier is <=
   the correct multiplier. Hence, lambda will be increased.  The slope of
   the dual function, neglecting b, starts out larger than b. We stop
   when we reach b. We assume that a >= 0, so that as lambda increases,
   x_i (lambda) decreases. Hence, the only bound variables that can become
   free are those with x_i (lambda) >= 1 */

#include "Napsack.hpp"

namespace SuiteSparse_Mongoose
{

Double QPnapup      /* return lambda */
(
    Double *x,              /* holds y on input, not modified */
    Int n,                  /* size of x */
    Double lambda,          /* initial guess for the shift */
    Double *a,              /* input constraint vector */
    Double b,               /* input constraint scalar */
    Double *breakpts,       /* break points */
    Int *bound_heap,        /* work array */
    Int *free_heap          /* work array */
)
{
    Int i, k, e, maxsteps, n_bound, n_free;
    Double ai, asum, a2sum, minbound, minfree, new_break, s, t, xi;

    minbound = INFINITY;
    minfree = INFINITY;

    /* -------------------------------------------------------------- */
    /* construct the heaps */
    /* -------------------------------------------------------------- */

    n_bound = 0;
    n_free = 0;
    asum = 0.;
    a2sum = 0.;
    if (a == NULL)
    {
        for (i = 0; i < n; i++)
        {
            xi = x[i] - lambda;
            if (xi > 1.)
            {
                n_bound++;
                bound_heap[n_bound] = i;
                asum++;
                t = x[i] - 1.;
                minbound = MONGOOSE_MIN2 (minbound, t);
                breakpts[i] = t;
            }
            else if (xi > 0.)
            {
                n_free++;
                free_heap[n_free] = i;
                asum += x[i];
                a2sum++;
                t = x[i];
                minfree = MONGOOSE_MIN2 (minfree, t);
                breakpts[i] = t;
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            ai = a[i];
            xi = x[i] - ai * lambda;
            if (xi > 1.)
            {
                n_bound++;
                bound_heap[n_bound] = i;
                asum += ai;
                t = (x[i] - 1.) / ai;
                minbound = MONGOOSE_MIN2 (minbound, t);
                breakpts[i] = t;
            }
            else if (xi > 0.)
            {
                n_free++;
                free_heap[n_free] = i;
                asum += x[i] * ai;
                a2sum += ai * ai;
                t = x[i] / ai;
                minfree = MONGOOSE_MIN2 (minfree, t);
                breakpts[i] = t;
            }
        }
    }

    maxsteps = 2 * n + 1;
    for (k = 1; k <= maxsteps; k++)
    {
        /*------------------------------------------------------------------- */
        /* check to see if zero slope achieved without changing the free set  */
        /* remember that the slope must always be adjusted by b               */
        /*------------------------------------------------------------------- */
        new_break = MONGOOSE_MIN2 (minfree, minbound);
        s = asum - new_break * a2sum;
        if ((s <= b) || (new_break == INFINITY)) /* done */
        {
            if (a2sum != 0.)
                lambda = (asum - b) / a2sum;
            return (lambda);
        }
        lambda = new_break;

        if (k == 1)
        {
            QPminheap_build(free_heap, n_free, breakpts);
            QPminheap_build(bound_heap, n_bound, breakpts);
        }

        /* -------------------------------------------------------------- */
        /* update the heaps */
        /* -------------------------------------------------------------- */

        if (n_free > 0)
        {
            if (a == NULL)
            {
                while (breakpts[e = free_heap[1]] <= lambda)
                {
                    a2sum--;
                    asum -= x[e];
                    n_free = QPminheap_delete(free_heap, n_free, breakpts);
                    if (n_free == 0)
                        break;
                }
            }
            else
            {
                while (breakpts[e = free_heap[1]] <= lambda)
                {
                    ai = a[e];
                    a2sum -= ai * ai;
                    asum -= ai * x[e];
                    n_free = QPminheap_delete(free_heap, n_free, breakpts);
                    if (n_free == 0)
                    {
                        a2sum = 0.;
                        break;
                    }
                }
            }
        }

        if (n_bound > 0)
        {
            if (a == NULL)
            {
                while (breakpts[e = bound_heap[1]] <= lambda)
                {
                    n_bound = QPminheap_delete(bound_heap, n_bound, breakpts);
                    a2sum++;
                    asum += x[e] - 1.;
                    breakpts[e] = x[e];
                    n_free = QPminheap_add(e, free_heap, breakpts, n_free);
                    if (n_bound == 0)
                        break;
                }
            }
            else
            {
                while (breakpts[e = bound_heap[1]] <= lambda)
                {
                    n_bound = QPminheap_delete(bound_heap, n_bound, breakpts);
                    ai = a[e];
                    a2sum += ai * ai;
                    asum += ai * (x[e] - 1.);
                    breakpts[e] = x[e] / ai;
                    n_free = QPminheap_add(e, free_heap, breakpts, n_free);
                    if (n_bound == 0)
                        break;
                }
            }
        }

        /*------------------------------------------------------------------- */
        /* get the biggest entry in each heap */
        /*------------------------------------------------------------------- */
        minfree = (n_free > 0 ? breakpts[free_heap[1]] : INFINITY);
        minbound = (n_bound > 0 ? breakpts[bound_heap[1]] : INFINITY);
    }

    /* this should not happen */
    assert(false);
    lambda = MONGOOSE_ZERO;
    return lambda;
}

}