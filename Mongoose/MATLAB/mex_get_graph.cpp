#include "gp_mex.hpp"

namespace SuiteSparse_Mongoose
{

Graph *mex_get_graph
(
    const mxArray *Gmatlab, /* The sparse matrix            */
    const mxArray *Amatlab  /* The real-valued node weights */
)
{
    Graph *returner = (Graph*) calloc(1, sizeof(Graph));
    if(!returner) return NULL;
    new (returner) Graph();

    Int n = returner->n = mxGetN(Gmatlab);
    Int *Gp = returner->p = (Int*) mxGetJc(Gmatlab);
    Int *Gi = returner->i = (Int*) mxGetIr(Gmatlab);
    Weight *Gx = returner->x = (Weight*) mxGetPr(Gmatlab);
    Int nz = returner->nz = Gp[n];

    /* Read node weights from matlab into the problem. */
    if(Amatlab != NULL)
    {
        returner->w = (Weight*) mxGetPr(Amatlab);
    }
    else
    {
        returner->w = (Weight*) malloc(n * sizeof(Weight));
        for(Int k=0; k<n; k++) returner->w[k] = 1.0;
    }

    return returner;
}

}