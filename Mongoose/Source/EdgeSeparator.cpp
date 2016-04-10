
#include "Mongoose_internal.hpp"
#include "Coarsening.hpp"
#include "GuessCut.hpp"
#include "Refining.hpp"
#include "Waterdance.hpp"

namespace SuiteSparse_Mongoose
{

bool initialize(Graph *G, Options *O);

/* The input must be a single connected component. */
void ComputeEdgeSeparator(Graph *G, Options *O)
{
    /* Check inputs */
    if(!G || !O) return;

    /* Finish initialization */
    if(!initialize(G, O)) return;

    /* Keep track of what the current graph is at any stage */
    Graph *current = G;

    /* If we need to coarsen the graph, do the coarsening. */
    while(current->n >= O->coarsenLimit)
    {
        match(current, O);
        Graph *next = coarsen(current, O);

        /* If we ran out of memory during coarsening, unwind the stack. */
        if(!next)
        {
             while(current != G)
             {
                  next = current->parent;
                  current->~Graph();
                  SuiteSparse_free(current);
                  current = next;
             }
             return;
        }

        current = next;
    }

    /* 
     * Generate a guess cut and do FM refinement.
     * On failure, unwind the stack.
     */
    if(!guessCut(current, O))
    {
         while(current != G)
         {
              Graph *next = current->parent;
              current->~Graph();
              SuiteSparse_free(current);
              current = next;
         }
         return;
    }

    /*
     * Refine the guess cut back to the beginning.
     */
    while(current->parent != NULL)
    {
        current = refine(current, O);
        waterdance(current, O);
    }
}

/* Finish the initialization of the top level graph. */
bool initialize(Graph *G, Options *O)
{
    Int n = G->n;
    Int *Gp = G->p;
    Weight *Gx = G->x;
    Weight *Gw = G->w;

    G->cn = 0;
    G->matching = (Int*) SuiteSparse_calloc(n, sizeof(Int));
    G->matchmap = (Int*) SuiteSparse_calloc(n, sizeof(Int));
    G->invmatchmap = (Int*) SuiteSparse_malloc(n, sizeof(Int));
    G->matchtype = (Int*) SuiteSparse_calloc(n, sizeof(Int));
    G->mark = (Int*) SuiteSparse_calloc(n, sizeof(Int));
    G->markValue = 1;

    G->partition = (bool*) SuiteSparse_malloc(n, sizeof(bool));
    G->bhIndex = (Int*) SuiteSparse_calloc(n, sizeof(Int));
    G->bhHeap[0] = (Int*) SuiteSparse_malloc(n, sizeof(Int));
    G->bhHeap[1] = (Int*) SuiteSparse_malloc(n, sizeof(Int));
    G->vertexGains = (Weight*) SuiteSparse_malloc(G->n, sizeof(Weight));
    G->externalDegree = (Int*) SuiteSparse_calloc(n, sizeof(Int));

    /* Check memory and abort if necessary. */
    if (!G->matching || !G->matchmap || !G->invmatchmap || !G->matchtype || 
        !G->mark || !G->partition || !G->bhIndex || !G->bhHeap[0] || !G->bhHeap[1] ||
        !G->vertexGains || !G->externalDegree) 
    {
        G->matching = (Int*) SuiteSparse_free(G->matching);
        G->matchmap = (Int*) SuiteSparse_free(G->matchmap);
        G->invmatchmap = (Int*) SuiteSparse_free(G->invmatchmap);
        G->matchtype = (Int*) SuiteSparse_free(G->matchtype);
        G->mark = (Int*) SuiteSparse_free(G->mark);
        G->partition = (bool*) SuiteSparse_free(G->partition);
        G->bhIndex = (Int*) SuiteSparse_free(G->bhIndex);
        G->bhHeap[0] = (Int*) SuiteSparse_free(G->bhHeap[0]);
        G->bhHeap[1] = (Int*) SuiteSparse_free(G->bhHeap[1]);
        G->vertexGains = (Weight*) SuiteSparse_free(G->vertexGains);
        G->externalDegree = (Int*) SuiteSparse_free(G->externalDegree);
        return false;
    }

    /* Compute worst-case gains, and compute X. */
    Weight X = 0.0, W = 0.0;
    Weight *gains = G->vertexGains;
    for(Int k=0; k<n; k++)
    {
        W += Gw[k];
        Weight sumEdgeWeights = 0.0;

        for(Int p=Gp[k]; p<Gp[k+1]; p++) sumEdgeWeights += Gx[p];

        gains[k] = -sumEdgeWeights;
        X += sumEdgeWeights;
    }
    G->X = X;
    G->W = W;
    G->H = 2.0 * X;
    return true;
}

}