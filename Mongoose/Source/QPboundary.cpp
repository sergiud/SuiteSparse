/* ========================================================================== */
/* === QPboundary =========================================================== */
/* ========================================================================== */

/*
Move all components of x to boundary of the feasible region

    0 <= x <= 1, a'x = b, lo <= b <= hi

while decreasing the cost function. The algorithm has the following parts

1. For each i in the free set, see if x_i can be feasibly pushed to either
   boundary while decreasing the cost.

2. For each i in the bound set, see if x_i can be feasibly flipped to opposite
   boundary while decreasing the cost.

3. For each i in the free list with a_{ij} = 0 and with j free,
   move either x_i or x_j to the boundary while decreasing
   the cost. The adjustments has the form x_i = s/a_i and x_j = -s/a_j 
   where s is a scalar factor. These adjustments must decrease cost.

4. For the remaining i in the free list, take pair x_i and x_j and
   apply adjustments of the same form as in #2 above to push at least one
   component to boundary. The quadratic terms can only decrease the
   cost function. We choose the sign of s such that g_i x_i + g_j x_j <= 0.
   Hence, this adjustment cannot increase the cost.
*/

/* ========================================================================== */

#include "GradProj.hpp"

namespace SuiteSparse_Mongoose
{

void QPboundary
(
    Graph *G,
    Options *O,
    QPDelta *QP
)
{
    /* ---------------------------------------------------------------------- */
    /* Step 0. read in the needed arrays                                      */
    /* ---------------------------------------------------------------------- */

    /* input and output */
    Int nf = QP->nf ;
    if(nf == 0) return;

    Double pert0 = 0.5 - 1e-10; //QP->Parm->boundary_pert0 ;
    Double pert1 = 0.5 + 1e-10; //QP->Parm->boundary_pert1 ;
    Double *x = QP->x ;         /* current estimate of solution */
    Double *grad = QP->g ;      /* gradient at current x */
    Int ib = QP->ib ;           /* ib = +1, -1, or 0 if b = hi, lo, or lo < b < hi */
    Int *ix = QP->ix ;          /* ix_i = +1, -1, or 0 if x_i = 1, 0, or 0 < x_i < 1*/
    Int *LinkUp = QP->LinkUp ;  /* linked list for free indices */
    Int *LinkDn = QP->LinkDn ;  /* linked list, LinkDn [LinkUp [i] = i*/
    Double b = QP->b ;          /* current value for a'x */

    /* problem specification */
    Int n  = G->n ;             /* problem dimension */
    Double *Ex = G->x ;         /* numerical values for edge weights */
    Int *Ei = G->i ;            /* adjacent vertices for each node */
    Int *Ep = G->p ;            /* points into Ex or Ei */
    Double *a  = G->w ;         /* a'x = b, lo <= b <= hi */
    Double lo = G->W * (O->targetSplit <= 0.5 ? O->targetSplit : 1 - O->targetSplit);
    Double hi = G->W * (O->targetSplit >= 0.5 ? O->targetSplit : 1 - O->targetSplit);
    Int *mark = G->mark;
    Int markValue = G->markValue;

    /* work array */
    Double *D  = QP->D ;   /* diagonal of quadratic */
//    Int *w0 = QP->wi0 ; /* Int size n, always zero */

    /* ---------------------------------------------------------------------- */
    /* Step 1. if lo < b < hi, then for each free j,                          */
    /*         see if x_j can be pushed to 0 or 1                             */
    /* ---------------------------------------------------------------------- */

    for(Int k=LinkUp[n]; (k<n) && (ib==0); k=LinkUp[k])
    {
        Double s, ak = a[k];
        if (grad[k] > 0.0) /* decrease x_j */
        {
            s = (b - lo) / ak;
            if(s < x[k])  /* lo is reached */
            {
                ib = -1;
                b = lo;
                x[k] -= s;
            }
            else          /* lo not reached */
            {
                s = x[k];
                x[k] = 0.;
                ix[k] = -1;
            }
        }
        else /* increase x_j */
        {
            s = (b - hi) / ak;
            if (s < x[k] - 1.) /* hi not reached */
            {
                s = x[k] - 1.;
                x[k] = 1.;
                ix[k] = +1;
            }
            else /* hi is reached */
            {
                ib = +1;
                b = hi;
                x[k] -= s;
            }
        }

        if(ib == 0) /* x_j hits boundary */
        {
            nf--;
            Int h = LinkUp[k];
            Int g = LinkDn[k];
            LinkUp[g] = h;
            LinkDn[h] = g;
            b -= s * ak;
        }

        for(Int p = Ep[k]; p < Ep[k+1]; p++)
        {
            grad[Ei[p]] += s * Ex[p];
        }
        grad[k] += s * D[k];
    }

    /* ---------------------------------------------------------------------- */
    /* Step 2. Examine flips of x_j from 0 to 1 or from 1 to 0 */
    /* ---------------------------------------------------------------------- */

    for(Int k=0; k<n; k++)
    {
        Int ixj = ix[k];
        if (ixj == 0) continue;

        Double ak = a[k];
        if (ixj > 0) /* try changing x_j from 1 to 0 */
        {
            if (b - ak >= lo)
            {
                if (0.5 * D[k] + grad[k] >= 0) /* flip lowers cost */
                {
                    b -= ak;
                    ib = (b == lo ? -1 : 0);
                    x[k] = MONGOOSE_ZERO;
                    ix[k] = -1;
                }
            }
        }
        else /* try changing x_j from 0 to 1 */
        {
            if (b + ak <= hi)
            {
                if (grad[k] - 0.5 * D[k] <= 0) /* flip lowers cost */
                {
                    b += ak;
                    ib = (b == hi ? 1 : 0);
                    x[k] = MONGOOSE_ONE;
                    ix[k] = +1;
                }
            }
        }

        if (ixj != ix[k])
        {
            if (ixj == 1) /* x [j] was 1, now it is 0 */
            {
                for(Int p = Ep[k]; p < Ep[k+1]; p++)
                {
                    grad[Ei[p]] += Ex[p];
                }
                grad[k] += D[k];
            }
            else /* x [j] was 0, now it is 1 */
            {
                for(Int p = Ep[k]; p < Ep[k+1]; p++)
                {
                    grad[Ei[p]] -= Ex[p];
                }
                grad[k] -= D[k];
            }
        }
    }

    if (nf == 0)
    {
        QP->nf = nf;
        QP->b = b;
        QP->ib = ib;
        return;
    }

    /* ---------------------------------------------------------------------- */
    /* Step 3. Search for a_{ij} = 0 in the free index set */
    /* ---------------------------------------------------------------------- */

    for(Int j = LinkUp[n]; j<n; j = LinkUp[j])
    {
        Int m = 1;
        for(Int p = Ep[j]; p < Ep[j+1]; p++)
        {
            if(ix[Ei[p]] == 0) m++;
        }
        if (m == nf) continue;

        /* ---------------------------------------------------------------------- */
        /* otherwise there exist i and j free with a_{ij} = 0, scatter Ei */
        /* ---------------------------------------------------------------------- */

        for(Int p = Ep[j]; p < Ep[j+1]; p++)
        {
            MONGOOSE_MARK(Ei[p]);
//            w0[Ei[p]] = 1;
        }
        MONGOOSE_MARK(j);
//        w0[j] = 1;

        for (Int i = LinkUp[n]; i < n; i = LinkUp[i])
        {
//            if (w0[i] == 0) /* a_{ij} = 0 */
            if(!MONGOOSE_MARKED(i))
            {
                Double aj = a[j];
                Double ai = a[i];
                Double xi = x[i];
                Double xj = x[j];

                /* cost change if x_j increases dx_j = s/a_j, dx_i = s/a_i */
                Double s;
                Int bind1, bind2;
                if (aj * (MONGOOSE_ONE - xj) < ai * xi) /* x_j hits upper bound */
                {
                    s = aj * (MONGOOSE_ONE - xj);
                    bind1 = 1;
                }
                else /* x_i hits lower bound */
                {
                    s = ai * xi;
                    bind1 = 0;
                }
                Double dxj = s / aj;
                Double dxi = -s / ai;
                Double c1 = (grad[j] - .5 * D[j] * dxj) * dxj + (grad[i] - .5 * D[i] * dxi) * dxi;

                /* cost change if x_j decreases dx_j = s/a_j, dx_i = s/a_i */
                if (aj * xj < ai * (MONGOOSE_ONE - xi)) /* x_j hits lower bound */
                {
                    s = -aj * xj;
                    bind2 = -1;
                }
                else /* x_i hits upper bound */
                {
                    s = -ai * (MONGOOSE_ONE - xi);
                    bind2 = 0;
                }
                dxj = s / aj;
                dxi = -s / ai;
                Double c2 = (grad[j] - 0.5 * D[j] * dxj) * dxj + (grad[i] - 0.5 * D[i] * dxi) * dxi;
                if (c1 < c2) /* increase x_j */
                {
                    if (bind1 == 1) /* x_j = 1 */
                    {
                        dxj = 1. - xj;
                        dxi = -aj * dxj / ai;
                        x[j] = 1.;
                        x[i] += dxi;
                        ix[j] = +1; /* j is bound at 1 */
                    }
                    else /* x_i = 0 */
                    {
                        dxi = -xi;
                        dxj = -ai * dxi / aj;
                        x[i] = 0.;
                        x[j] += dxj;
                        ix[i] = -1; /* i is bound at 0 */
                    }
                }
                else
                {
                    if (bind2 == -1) /* x_j = 0 */
                    {
                        bind1 = 1; /* j is bound, not i */
                        x[j] = 0.;
                        x[i] += dxi;
                        ix[j] = -1; /* j is bound at 0 */
                    }
                    else /* x_i = 1 */
                    {
                        bind1 = 0; /* j not bound */
                        x[i] = 1;
                        x[j] += dxj;
                        ix[i] = +1; /* i is bound at 0 */
                    }
                }

                for (Int p = Ep[j]; p < Ep[j+1]; p++) grad[Ei[p]] -= Ex[p] * dxj;
                for (Int p = Ep[i]; p < Ep[i+1]; p++) grad[Ei[p]] -= Ex[p] * dxi;
                grad[j] -= D[j] * dxj;
                grad[i] -= D[i] * dxi;
                nf--;

                if (bind1) /* remove j from free list */
                {
                    Int h = LinkUp[j];
                    Int g = LinkDn[j];
                    LinkUp[g] = h;
                    LinkDn[h] = g;
                    break;
                }
                else /* remove i from free list */
                {
                    Int h = LinkUp[i];
                    Int g = LinkDn[i];
                    LinkUp[g] = h;
                    LinkDn[h] = g;
                    continue;
                }
            }
        }

        markValue++;
#if 0
        /* ---------------------------------------------------------------------- */
        /* restore w0 to zero */
        /* ---------------------------------------------------------------------- */
        for(Int k = Ep[j]; k<Ep[j+1]; k++) w0[Ei[k]] = 0;
        w0[j] = 0;
#endif
    }

    /* ---------------------------------------------------------------------- */
    /* Step 4. dxj = s/aj, dxi = -s/ai, choose s with g_j dxj + g_i dxi <= 0 */
    /* ---------------------------------------------------------------------- */

    /* free variables:0 < x_j < 1 */

    Int j;
    for(j = LinkUp[n]; j < n; j = LinkUp[j]) /* free variables:0 < x_j < 1 */
    {
        /* choose s so that first derivative terms decrease */

        Int i = LinkUp[j];
        if(i == n) break;
        Double ai = a[i];
        Double aj = a[j];
        Double xi = x[i];
        Double xj = x[j];

        Int bind1;
        Double dxj, dxi, s = grad[j] / aj - grad[i] / ai;
        if (s < MONGOOSE_ZERO) /* increase x_j */
        {
            if (aj * (MONGOOSE_ONE - xj) < ai * xi) /* x_j hits upper bound */
            {
                dxj = MONGOOSE_ONE - xj;
                dxi = -aj * dxj / ai;
                x[j] = MONGOOSE_ONE;
                x[i] += dxi;
                ix[j] = +1;
                bind1 = 1; /* x_j is bound */
            }
            else /* x_i hits lower bound */
            {
                dxi = -xi;
                dxj = -ai * dxi / aj;
                x[i] = 0.;
                x[j] += dxj;
                ix[i] = -1;
                bind1 = 0; /* x_i is bound */
            }
        }
        else /* decrease x_j */
        {
            if (aj * xj < ai * (1. - xi)) /* x_j hits lower bound */
            {
                dxj = -xj;
                dxi = -aj * dxj / ai;
                x[j] = 0;
                x[i] += dxi;
                ix[j] = -1;
                bind1 = 1; /* x_j is bound */
            }
            else /* x_i hits upper bound */
            {
                dxi = 1 - xi;
                dxj = -ai * dxi / aj;
                x[i] = 1;
                x[j] += dxj;
                ix[i] = +1;
                bind1 = 0; /* x_i is bound */
            }
        }

        for(Int k = Ep[j]; k < Ep[j+1]; k++) grad[Ei[k]] -= Ex[k] * dxj;
        for(Int k = Ep[i]; k < Ep[i+1]; k++) grad[Ei[k]] -= Ex[k] * dxi;
        grad[j] -= D[j] * dxj;
        grad[i] -= D[i] * dxi;

        if (bind1) /* j is bound */
        {
            Int h = LinkUp[j];
            Int g = LinkDn[j];
            LinkUp[g] = h;
            LinkDn[h] = g;
        }
        else /* i is bound */
        {
            Int h = LinkUp[i];
            Int g = LinkDn[i];
            LinkUp[g] = h;
            LinkDn[h] = g;
            j = LinkDn[j]; /* j is still free, repeat it */
        }
        nf--;
    }

    if (nf == 1) /* j is free, optimize over x [j] */
    {
        Int bind1 = 0;
        Double aj = a[j];
        Double dxj = (hi - b) / aj;
        if (dxj < MONGOOSE_ONE - x[j])
        {
            bind1 = 1;
        }
        else
        {
            dxj = MONGOOSE_ONE - x[j];
        }

        Int bind2 = 0;
        Double dxi = (lo - b) / aj;
        if (dxi > -x[j])
        {
            bind2 = 1;
        }
        else
        {
            dxi = -x[j];
        }

        Double c1 = (grad[j] - 0.5 * D[j] * dxj) * dxj;
        Double c2 = (grad[j] - 0.5 * D[j] * dxi) * dxi;
        if (c1 <= c2) /* x [j] += dxj */
        {
            if (bind1)
            {
                x[j] += dxj;
                ib = +1;
                b = hi;
            }
            else
            {
                nf--;
                LinkUp[n] = n;
                LinkDn[n] = n;
                x[j] = 1.;
                ix[j] = 1;
                b += dxj * aj;
            }
        }
        else /* x [j] += dxi */
        {
            dxj = dxi;
            if (bind2)
            {
                x[j] += dxj;
                ib = -1;
                b = lo;
            }
            else
            {
                nf--;
                LinkUp[n] = n;
                LinkDn[n] = n;
                x[j] = 0.;
                ix[j] = -1;
                b += dxj * aj;
            }
        }

        if (dxj != MONGOOSE_ZERO)
        {
            for(Int p = Ep[j]; p < Ep[j+1]; p++)
            {
                grad[Ei[p]] -= Ex[p] * dxj;
            }
            grad[j] -= D[j] * dxj;
        }

        if((x[j] >= pert1) || ((aj == MONGOOSE_ONE) && (x[j] > 0.5)))
        {
            nf--;
            LinkUp[n] = n;
            LinkDn[n] = n;
            x[j] = 1.;
            ix[j] = +1;
        }
        else if ((x[j] <= pert0) || ((aj == 1.) && (x[j] < .5)))
        {
            nf--;
            LinkUp[n] = n;
            LinkDn[n] = n;
            x[j] = 0.;
            ix[j] = -1;
        }
    }
    QP->nf = nf;
    QP->b = b;
    QP->ib = ib;

    G->markValue = markValue + 1;
}

}