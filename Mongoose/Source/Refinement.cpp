
#include "Refining.hpp"
#include "Matching.hpp"
#include "BoundaryHeap.hpp"
#include "FidduciaMattheyes.hpp"

namespace SuiteSparse_Mongoose
{

Graph *refine(Graph *G, Options *O)
{
    Graph *P = G->parent;

    Int cn = G->n;
//    Int *Gp = G->p;
//    Int *Gi = G->i;
    bool *cPartition = G->partition;
//    Int *bhIndex = G->bhIndex;

//    Int fn = P->n;
//    Int *Pp = P->p;
//    Int *Pi = P->i;
    Int *invmatchmap = P->invmatchmap;
    Int *matching = P->matching;
    bool *fPartition = P->partition;
    Weight *fGains = P->vertexGains;
    Int *fExternalDegree = P->externalDegree;

    /* Transfer cut costs and partition details upwards. */
    P->heuCost = G->heuCost;
    P->cutCost = G->cutCost;
    P->W0 = G->W0;
    P->W1 = G->W1;
    P->imbalance = G->imbalance;

    /* For each vertex in the coarse graph. */
    for(Int k=0; k<cn; k++)
    {
        /* Load up the inverse matching */
        Int v[3] = {-1, -1, -1};
        v[0] = invmatchmap[k];
        v[1] = MONGOOSE_GETMATCH(v[0]);
        if(v[0] == v[1]){ v[1] = -1; }
        else
        {
            v[2] = MONGOOSE_GETMATCH(v[1]);
            if(v[0] == v[2]){ v[2] = -1; }
        }

        /* Transfer the partition choices to the fine level. */
        bool cp = cPartition[k];
        for(Int i=0; i<3 && v[i] != -1; i++)
        {
            Int vertex = v[i];
            fPartition[vertex] = cp;
        }
    }

    /* See if we can relax the boundary constraint and recompute gains for
     * vertices on the boundary.
     * NOTE: For this, we only need to go through the set of vertices that
     * were on the boundary in the coarse representation. */
    for(Int h=0; h<2; h++)
    {
        /* Get the appropriate heap's data. */
        Int *heap = G->bhHeap[h];
        Int size = G->bhSize[h];

        /* Go through all the boundary nodes. */
        for(Int hpos=0; hpos<size; hpos++)
        {
            /* Get the coarse vertex from the heap. */
            Int k = heap[hpos];

            /* Load up the inverse matching */
            Int v[3] = {-1, -1, -1};
            v[0] = invmatchmap[k];
            v[1] = MONGOOSE_GETMATCH(v[0]);
            if(v[0] == v[1]){ v[1] = -1; }
            else
            {
                v[2] = MONGOOSE_GETMATCH(v[1]);
                if(v[0] == v[2]){ v[2] = -1; }
            }

            /* Relax the boundary constraint. */
            for(Int i=0; i<3 && v[i] != -1; i++)
            {
                Int vertex = v[i];

                Weight gain;
                Int externalDegree;
                calculateGain(P, O, vertex, &gain, &externalDegree);

                /* Only add relevant vertices to the boundary heap. */
                if(externalDegree > 0)
                {
                    fExternalDegree[vertex] = externalDegree;
                    fGains[vertex] = gain;
                    bhInsert(P, vertex);
                }
            }
        }
    }

    /* Now that we're done with the coarse graph, we can release it. */
    G->~Graph();
    G = (Graph*) SuiteSparse_free(G);

// writeDot(P, O, "refined", CutSet);

    return P;
}

}
