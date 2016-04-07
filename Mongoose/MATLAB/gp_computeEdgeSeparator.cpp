
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
    int i, gtype, values;

    const char* usage = "Usage: partition = gp_computeEdgeSeparator(G, (O, A))";
    if(nargout != 1 || nargin < 1 || nargin > 3) mexErrMsgTxt(usage);

    const mxArray *matGraph = pargin[0];
    const mxArray *matOptions = (nargin >= 2 ? pargin[1] : NULL);
    const mxArray *matNodeWeights = (nargin >= 3 ? pargin[2] : NULL);

    /* Get the graph from the matlab inputs. */
    Graph *G = mex_get_graph(matGraph, matNodeWeights);
    if(!G)
        mexErrMsgTxt("Unable to get Graph struct");

    /* Get the options from the matlab inputs. */
    Options *O = mex_get_options(matOptions);
    if(!O)
        mexErrMsgTxt("Unable to get Options struct");

    ComputeEdgeSeparator(G, O);

#if 0
    /* Compute the edge separator */
    writeDot(G, O, "initial");
    computeEdgeSeparator(G, O);
    writeDot(G, O, "final", CutSet);
#endif

    /* Copy the partition choices back to matlab. */
    pargout[0] = gp_mex_put_logical(G->partition, G->n) ;

    G->p = NULL;
    G->i = NULL;
    G->x = NULL;
    if(matNodeWeights != NULL) G->w = NULL;

    /* Cleanup */
    G->~Graph();
    free(G);
    O->~Options();
    free(O);
}
