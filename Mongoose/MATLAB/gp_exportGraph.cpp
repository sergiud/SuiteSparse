
#include "gp_mex.hpp"

using namespace SuiteSparse_Mongoose;

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    cs Amatrix;

    const char* usage = "Usage: gp_export(A)";
    if(nargout != 0 || nargin != 1) mexErrMsgTxt(usage);

    Graph *G = mex_get_graph(pargin[0]);
    if(!G) return;
    Int n = G->n;
    Int *Gp = G->p;
    Int *Gi = G->i;
    Weight *Gx = G->x;

    FILE *troll = fopen("troll.mtx", "w");
    fprintf(troll, "%ld %ld %ld\n", n, n, Gp[n]);
    for(Int k=0; k<n; k++)
    {
        for(Int p=Gp[k]; p<Gp[k+1]; p++)
        {
            fprintf(troll, "%ld %ld %f\n", k+1, Gi[p]+1, Gx[p]);
        }
    }
    fclose(troll);

    /* Remove linked data structures before deleting the graph. */
    G->p = NULL;
    G->i = NULL;
    G->x = NULL;
    G->~Graph();
    free(G);
}

