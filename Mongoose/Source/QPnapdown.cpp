/* ========================================================================== */
/* === QPnapdown ============================================================ */
/* ========================================================================== */

/* Find x that minimizes ||x-y|| while satisfying the constraints
   0 <= x <= 1, a'x = b. If a = NULL, then we assume that a = 1.
   The algorithm is described in the napsack comments.
   It is assumed that the starting guess lambda for the dual multiplier is >=
   the correct multiplier. Hence, lambda will be decreased.  The slope of the
   dual function, neglecting b, starts out smaller than b. We stop
   when we reach b. We assume that a >= 0, so that as lambda decreases,
   x_i (lambda) increases. Hence, the only bound variables that can become
   free are those with x_i (lambda) <= 0 */

#include "Napsack.hpp"

namespace SuiteSparse_Mongoose
{

Double QPnapdown            /* return lambda */
(
    Double *x,	            /* holds y on input, not modified */
    Int n,	            /* size of x */
    Double lambda,          /* initial guess for the shift */
    Double *a,              /* input constraint vector */
    Double b,        	    /* input constraint scalar */
    Double *breakpts,       /* break points */
    Int *bound_heap,        /* work array */
    Int *free_heap          /* work array */
)
{
    Int i, k, e, maxsteps, n_bound, n_free;
    Double ai, asum, a2sum, maxbound, maxfree, new_break, s, t, xi;

    maxbound = -INFINITY;
    maxfree = -INFINITY;

    /* -------------------------------------------------------------- */
    /* construct the heaps */
    /* -------------------------------------------------------------- */

    n_bound = 0;
    n_free = 0;
    asum = MONGOOSE_ZERO;
    a2sum = MONGOOSE_ZERO;
    if (a == NULL)
    {
        for (i = 0; i < n; i++)
        {
            xi = x[i] - lambda;
            if (xi < MONGOOSE_ZERO)
            {
                n_bound++;
                bound_heap[n_bound] = i;
                t = x[i];
                maxbound = MONGOOSE_MAX2(maxbound, t);
                breakpts[i] = t;
            }
            else if (xi < MONGOOSE_ONE)
            {
                n_free++;
                free_heap[n_free] = i;
                t = x[i] - MONGOOSE_ONE;
                asum += x[i];
                a2sum++;
                maxfree = MONGOOSE_MAX2(maxfree, t);
                breakpts[i] = t;
            }
            else
                asum++;
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            ai = a[i];
            xi = x[i] - ai * lambda;
            if (xi < MONGOOSE_ZERO)
            {
                n_bound++;
                bound_heap[n_bound] = i;
                t = x[i] / ai;
                maxbound = MONGOOSE_MAX2(maxbound, t);
                breakpts[i] = t;
            }
            else if (xi < MONGOOSE_ONE)
            {
                n_free++;
                free_heap[n_free] = i;
                t = (x[i] - MONGOOSE_ONE) / ai;
                asum += x[i] * ai;
                a2sum += ai * ai;
                maxfree = MONGOOSE_MAX2(maxfree, t);
                breakpts[i] = t;
            }
            else
                asum += ai;
        }
    }

    /*------------------------------------------------------------------- */
    /* check to see if zero slope achieved without changing the free set  */
    /* remember that the slope must always be adjusted by b               */
    /*------------------------------------------------------------------- */

    maxsteps = 2 * n + 1;
    for (k = 1; k <= maxsteps; k++)
    {
        new_break = MONGOOSE_MAX2(maxfree, maxbound);
        s = asum - new_break * a2sum;
        if ((s >= b) || (new_break == -INFINITY)) /* done */
        {
            if(a2sum != MONGOOSE_ZERO) lambda = (asum - b) / a2sum;
            return lambda;
        }
        lambda = new_break;

        if (k == 1)
        {
            QPmaxheap_build(free_heap, n_free, breakpts);
            QPmaxheap_build(bound_heap, n_bound, breakpts);
        }

        /* -------------------------------------------------------------- */
        /* update the heaps */
        /* -------------------------------------------------------------- */

        if (n_free > 0)
        {
            if (a == NULL)
            {
                while (breakpts[e = free_heap[1]] >= lambda)
                {
                    a2sum--;
                    asum = asum + (MONGOOSE_ONE - x[e]);

                    n_free = QPmaxheap_delete(free_heap, n_free, breakpts);
                    if (n_free == 0)
                        break;
                }
            }
            else
            {
                while (breakpts[e = free_heap[1]] >= lambda)
                {
                    ai = a[e];
                    a2sum -= ai * ai;
                    asum = asum + ai * (MONGOOSE_ONE - x[e]);

                    n_free = QPmaxheap_delete(free_heap, n_free, breakpts);
                    if (n_free == 0)
                    {
                        a2sum = MONGOOSE_ZERO;
                        break;
                    }
                }
            }
        }

        if (n_bound > 0)
        {
            if (a == NULL)
            {
                while (breakpts[e = bound_heap[1]] >= lambda)
                {
                    n_bound = QPmaxheap_delete(bound_heap, n_bound, breakpts);
                    a2sum++;
                    asum = asum + x[e];
                    t = x[e] - MONGOOSE_ONE;
                    breakpts[e] = t;
                    n_free = QPmaxheap_add(e, free_heap, breakpts, n_free);
                    if (n_bound == 0)
                        break;
                }
            }
            else
            {
                while (breakpts[e = bound_heap[1]] >= lambda)
                {
                    n_bound = QPmaxheap_delete(bound_heap, n_bound, breakpts);
                    ai = a[e];
                    a2sum += ai * ai;
                    asum = asum + ai * x[e];
                    t = (x[e] - MONGOOSE_ONE) / ai;
                    breakpts[e] = t;
                    n_free = QPmaxheap_add(e, free_heap, breakpts, n_free);
                    if (n_bound == 0)
                        break;
                }
            }
        }

        /*------------------------------------------------------------------- */
        /* get the biggest entry in each heap */
        /*------------------------------------------------------------------- */

        maxfree = (n_free > 0 ? breakpts[free_heap[1]] : -INFINITY);
        maxbound = (n_bound > 0 ? breakpts[bound_heap[1]] : -INFINITY);
    }

    /* This should never happen */
    assert(false);
    lambda = MONGOOSE_ZERO;
    return lambda;
}

}