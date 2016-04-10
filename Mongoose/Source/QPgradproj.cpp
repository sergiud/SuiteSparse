/* ========================================================================== */
/* === QPgradproj =========================================================== */
/* ==========================================================================
 Apply gradient projection algorithm to the quadratic program which
 arises in graph partitioning:

 min (1-x)'(D+A)x subject to lo <= b <= hi, a'x = b, 0 <= x <= 1

 The gradient at the current point is provided as input, and the
 gradient is updated in each iteration.
 ============================================================================ */

#include "GradProj.hpp"
#include "Napsack.hpp"

namespace SuiteSparse_Mongoose
{

Double QPgradproj
(
    Graph *G,
    Options *O,
    QPDelta *QP
)
{
    /* ---------------------------------------------------------------------- */
    /* Unpack the relevant structures                                         */
    /* ---------------------------------------------------------------------- */
    Double tol = O->gradprojTol;
    Double *wx1 = QP->wx[0];    /* Double work array, used in napsack and main program     */
    Double *wx2 = QP->wx[1];    /* Double work array, used in napsack and as Dgrad in main */
    Double *wx3 = QP->wx[2];    /* Double work array, used in main for y - x               */
    Int *wi1 = QP->wi[0];       /* Integer work array, used in napsack and as C in main    */
    Int *wi2 = QP->wi[1];       /* Integer work array, used in napsack                     */

    /* Output and Input */
    Double *x = QP->x;             /* current estimate of solution            */
    Int nf = QP->nf;               /* number of i such that 0 < x_i < 1       */
    Int *ix = QP->ix;              /* ix_i = +1,-1, or 0 if x_i = 1,0, or 0 < x_i < 1 */
    Int *LinkUp = QP->LinkUp;      /* linked list for free indices            */
    Int *LinkDn = QP->LinkDn;      /* linked list, LinkDn [LinkUp [i]] = i    */
    Double *grad = QP->g;          /* gradient at current x                   */

    /* Unpack the problem's parameters. */
    Int n = G->n;                  /* problem dimension */
    Int *Ep = G->p;                /* points into Ex or Ei */
    Int *Ei = G->i;                /* adjacent vertices for each node */
    Weight *Ex = G->x;             /* edge weights */
    Weight *Ew = G->w;             /* node weights; a'x = b, lo <= b <= hi */
    Double lo = G->W * (O->targetSplit <= 0.5 ? O->targetSplit : 1 - O->targetSplit);
    Double hi = G->W * (O->targetSplit >= 0.5 ? O->targetSplit : 1 - O->targetSplit);

    Double *D = QP->D; /* diagonal of quadratic */

    /* gradient projection parameters */
    Int limit = O->gradprojIterationLimit; /* max number of iterations */

    /* work arrays */
    Double *y = wx1;
    Double *wx = wx2;
    Double *d = wx3;
    Double *Dgrad = wx;     /* gradient change       ; used in napsack as wx  */
    Int *C = wi1;           /* components of x change; used in napsack as wi1 */

    /* compute error, take step along projected gradient */
    Int ib = 0;             /* initialize ib so that lo < b < hi */
    Double lambda = MONGOOSE_ZERO;
    Int it = 0;
    Double err = INFINITY;

    /* TODO: ?? is this comment a correct interpretation ?? */
    /* Project x into feasible set. */
    for(Int k=0; k<n; k++)
    {
        ix[k] = (x[k] >= MONGOOSE_ONE ? 1 : x[k] <= MONGOOSE_ZERO ? -1 : 0);
    }

    while (err > tol)
    {
/* TODO: ?? is this comment a correct interpretation ?? */
        /* Moving in the gradient direction. */
        for(Int k=0; k<n; k++) y[k] = x[k] - grad[k];

        /* Run the napsack. */
        lambda = QPnapsack(y, n, lo, hi, Ew, lambda, ix, wx, wi1, wi2);

/* TODO: This can probably be done in the napsack to avoid O(n) loop. */
        /* Compute the maximum error. */
        err = -INFINITY;
        for(Int k=0; k<n; k++) err = MONGOOSE_MAX2(err, fabs(y[k]-x[k]));

        /* If we converged or got exhausted, save context and exit. */
        if ((err <= tol) || (it >= limit))
        {
            QP_GRADPROJ_saveContextAndExit();
        }

        it++;

        /* compute stepsize st = g_F'g_F/-g_F'(A+D)g_F */
        /* TODO: Can Dgrad be cleared after use to avoid O(n)? */
        for (Int k = 0; k < n; k++) Dgrad[k] = MONGOOSE_ZERO;
Int lastLink = -1;
        for (Int i = LinkUp[n]; i < n; i = LinkUp[i])
        {
/* INFINITE LOOP i=LinkUp[i] links to itself. */
if(i != lastLink)
{
    lastLink = i;
}
else
{
    abort();
}

            /* compute -(A+D)g_F */
            Double s = grad[i];
            for (Int p = Ep[i]; p < Ep[i+1]; p++)
            {
                Dgrad[Ei[p]] -= s * Ex[p];
            }
            Dgrad[i] -= s * D[i];
        }

        Double st_num = MONGOOSE_ZERO;
        Double st_den = MONGOOSE_ZERO;
        for (Int j = LinkUp[n]; j < n; j = LinkUp[j])
        {
            st_num += grad[j] * grad[j];
            st_den += grad[j] * Dgrad[j];
        }

        /* st = g_F'g_F/-g_F'(A+D)g_F unless the denominator <= 0 */
        if (st_den > MONGOOSE_ZERO)
        {
            Double st = MONGOOSE_MAX2 (st_num / st_den, 0.001);
            for (Int j = 0; j < n; j++) y[j] = x[j] - st * grad[j];
            lambda = QPnapsack(y, n, lo, hi, Ew, lambda, ix, wx, wi1, wi2);
        }

        /* otherwise st = 1 and y is as computed above */
        Int nc = 0; /* number of changes (number of j for which y_j != x_j) */
        Weight s = MONGOOSE_ZERO;
        for (Int j = 0; j < n; j++) Dgrad[j] = MONGOOSE_ZERO;
        for (Int j = 0; j < n; j++)
        {
            Double t = y[j] - x[j];
            if (t != MONGOOSE_ZERO)
            {
                d[j] = t;
                s += t * grad[j]; /* derivative in the direction y - x */
                C[nc] = j;
                nc++;
                for(Int p=Ep[j]; p<Ep[j+1]; p++)
                {
                    Dgrad[Ei[p]] -= Ex[p] * t;
                }
                Dgrad[j] -= D[j] * t;
            }
        }

        /* If directional derivative has wrong sign, save context and exit. */
        if(s >= MONGOOSE_ZERO)
        {
//printf("%d: Directional derivative has wrong sign (%f).\n", it, s);
            QP_GRADPROJ_saveContextAndExit();
        }

        Double t = MONGOOSE_ZERO;
        for (Int k = 0; k < nc; k++)
        {
            Int j = C[k];
            t += Dgrad[j] * d[j]; /* -dg'd */
        }

        if(s+t <= 0) /* min attained at y, slope at y <= 0 */
        {
            ib = (lambda > 0 ? 1 : lambda < 0 ? -1 : 0);
            for(Int k = 0; k < nc; k++)
            {
                Int j = C[k];
                Double yj = y[j];
                x[j] = yj;
                Int bind = -1; /* -1 = no change, 0 = free, +1 = bind */
                if(ix[j] > 0)
                {
                    if (yj == MONGOOSE_ZERO)
                    {
                        ix[j] = -1;
                    }
                    else
                    {
                        bind = 0; /* FREEJ */
                    }
                }
                else if (ix[j] < 0)
                {
                    if (yj == 1.0)
                    {
                        ix[j] = 1.0;
                    }
                    else
                    {
                        bind = 0;
                    }
                }
                else /* x_j currently free */
                {
                    if (yj == 1.0) /* x_j hits upper bound */
                    {
                        bind = 1;
                        ix[j] = 1.0;
                    }
                    else if (yj == MONGOOSE_ZERO) /* x_j hits lower bound */
                    {
                        bind = 1;
                        ix[j] = -1.0;
                    }
                    /* else x_j is free before and after step */
                }
                if(bind == 0)
                {
                    ix[j] = 0;
                    nf++;
                    Int m = LinkUp[n];
                    LinkUp[j] = m;
                    LinkUp[n] = j;
                    LinkDn[m] = j;
                    LinkDn[j] = n;
                }
                else if(bind == 1)
                {
                    nf--;
                    Int h = LinkUp[j];
                    Int g = LinkDn[j];
                    LinkUp[g] = h;
                    LinkDn[h] = g;
                }
            }
            for (Int j = 0; j < n; j++)
            {
                grad[j] += Dgrad[j];
            }
        }
        else /* partial step towards y, st < 1 */
        {
            if((ib > 0 && lambda <= 0) || (ib < 0 && lambda >= 0))
            {
                ib = 0;
            }

            Double st = -s / t;
            for (Int k = 0; k < nc; k++)
            {
                Int j = C[k];
                if(ix[j] != 0) /* x_j became free */
                {
                    ix[j] = 0;
                    nf++;
                    Int m = LinkUp[n];
                    LinkUp[j] = m;
                    LinkUp[n] = j;
                    LinkDn[m] = j;
                    LinkDn[j] = n;
                }

                /*  else x_j is free before and after step */
                x[j] += st * d[j];
            }

            for(Int k=0; k<n; k++) grad[k] += st * Dgrad[k];
        }
    }
    return err;
}

}