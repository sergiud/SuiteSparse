/* ========================================================================== */
/* === minheap ============================================================== */
/* ========================================================================== */

#include "Napsack.hpp"

namespace SuiteSparse_Mongoose
{

/* ========================================================================== */
/* === minheap_build ======================================================== */
/* ========================================================================== */

/* build a min heap in heap [1..nheap] */

void QPminheap_build
(
    Int *heap,	    /* on input, an unsorted set of elements */
    Int size,		/* number of elements to build into the heap */
    Double *x
)
{
    Int p;

    for (p = size / 2; p >= 1; p--)
    {
        QPminheapify(p, heap, size, x);
    }
}

/* ========================================================================== */
/* === minheap_delete ====================================================== */
/* ========================================================================== */

/* delete the top element in a min heap */

Int QPminheap_delete	/* return new size of heap */
(
    Int *heap,	/* containing indices into x, 1..n on input */
    Int size,		/* number of items in heap */
    Double *x	/* not modified */
)
{
    if (size <= 1)
    {
        return (0);
    }

    /* move element from the end of the heap to the top */
    heap[1] = heap[size];
    size--;
    QPminheapify(1, heap, size, x);
    return (size);
}


/* ========================================================================== */
/* === minheap_add ========================================================== */
/* ========================================================================== */
/* ========================================================================== */

/* add a new leaf to a min heap */

Int QPminheap_add
(
    Int leaf,   /* the new leaf */
    Int *heap,  /* size n, containing indices into x */
    Double *x,  /* not modified */
    Int nheap   /* number of elements in heap not counting new one */
)
{
    Int l, lnew, lold;
    Double xold, xnew;

    nheap++;
    lold = nheap;
    heap[lold] = leaf;
    xold = x[leaf];
    while (lold > 1)
    {
        lnew = lold / 2;
        l = heap[lnew];
        xnew = x[l];

        /* swap new and old */
        if (xnew > xold)
        {
            heap[lnew] = leaf;
            heap[lold] = l;
        }
        else
        {
            return (nheap);
        }

        lold = lnew;
    }
    return (nheap);
}

/* ========================================================================== */
/* === minheapify =========================================================== */
/* ========================================================================== */

/* heapify starting at node p.  On input, the heap at node p satisfies the */
/* heap property, except for heap [p] itself.  On output, the whole heap */
/* satisfies the heap property. */

void QPminheapify
(
    Int p,		/* start at node p in the heap */
    Int *heap,	/* size n, containing indices into x */
    Int size,		/* heap [ ... nheap] is in use */
    Double *x	/* not modified */
)
{
    Int left, right, e, hleft, hright;
    Double xe, xleft, xright;

    e = heap[p];
    xe = x[e];

    while(true)
    {
        left = p * 2;
        right = left + 1;

        if(right <= size)
        {
            hleft = heap[left];
            hright = heap[right];
            xleft = x[hleft];
            xright = x[hright];
            if (xleft < xright)
            {
                if (xe > xleft)
                {
                    heap[p] = hleft;
                    p = left;
                }
                else
                {
                    heap[p] = e;
                    return;
                }
            }
            else
            {
                if (xe > xright)
                {
                    heap[p] = hright;
                    p = right;
                }
                else
                {
                    heap[p] = e;
                    return;
                }
            }
        }
        else
        {
            if (left <= size)
            {
                hleft = heap[left];
                xleft = x[hleft];
                if (xe > xleft)
                {
                    heap[p] = hleft;
                    p = left;
                }
            }
            heap[p] = e;
            return;
        }
    }
}

/* ========================================================================== */
/* === minheap_check ======================================================== */
/* ========================================================================== */

/* ensure that heap [p..nheap] satisfies the heap property, and that all */
/* values of x are within the proper range. */

void QPminheap_check
(
    Int *heap,	/* vector of size n+1 */
    Double *x,	/* vector of size n */
    Int size,		/* # items in heap */
    Int n,
    Int p		/* start checking at heap [p] */
)
{
    Int *w, left, right, i, e, eleft, eright;
    Double xe;

    if (size == 0)
        return;

    /* allocate workspace */
    w = (Int*) malloc(n * sizeof(Int));

    for (i = 0; i < n; i++)
    {
        w[i] = 0;
    }

    for (i = p; i <= size; i++)
    {
        e = heap[i];

        w[e] = 1;

        left = i * 2;
        right = left + 1;

        /* get x [e] and make sure it's in range */
        xe = x[e];

        /* make sure that the node is <= left child */
        if (left <= size)
        {
            eleft = heap[left];
        }

        /* make sure that the node is <= right child */
        if (right <= size)
        {
            eright = heap[right];
        }
    }

    /* free workspace ] */
    free(w);
}

}