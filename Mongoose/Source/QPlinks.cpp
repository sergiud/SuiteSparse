
#include "GradProj.hpp"

namespace SuiteSparse_Mongoose
{

void QPlinks
(
    Graph *G,
    Options *O,
    QPDelta *QP          /* pointer to QPDelta structure  */
)
{
    /* Inputs */
    Double *x = QP->x;

    /* Unpack structures. */
    Int n = G->n;
    Int *Ep = G->p;
    Int *Ei = G->i;
    Weight *Ex = G->x;
    Weight *a = G->w;

    /* working array */
    Double *D = QP->D ;
    Int *ix = QP->ix ;
    Int *LinkUp = QP->LinkUp ;
    Int *LinkDn = QP->LinkDn ;
    Double *grad = QP->g ; /* gradient at current x */

    Int lastl = n;
    Int nf = 0;
    Double s = MONGOOSE_ZERO;

    for(Int k=0; k<n; k++) grad[k] = (0.5-x[k]) * D[k];
    for(Int k=0; k<n; k++)
    {
        Double xk = x[k];
        s += a[k] * xk;
        Weight r = 0.5 - xk;
        for(Int p=Ep[k]; p<Ep[k+1]; p++) grad[Ei[p]] += r * Ex[p];

        if(xk < MONGOOSE_ONE || xk > MONGOOSE_ZERO)
        {
            ix[k] = 0;
            LinkUp[lastl] = k;
            LinkDn[k] = lastl;
            lastl = k;
            nf++;
        }
        else
        {
            ix[k] = (xk >= MONGOOSE_ONE ? 1 : -1);
        }
    }

    LinkUp[lastl] = n;
    LinkDn[n] = lastl;
    QP->nf = nf;
    QP->b = s;

    Double lo = G->W * (O->targetSplit <= 0.5 ? O->targetSplit : 1 - O->targetSplit);
    Double hi = G->W * (O->targetSplit >= 0.5 ? O->targetSplit : 1 - O->targetSplit);
    QP->ib = (s <= lo ? -1 : s < hi ? 0 : 1);
}

}