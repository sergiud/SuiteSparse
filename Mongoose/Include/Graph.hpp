#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include "stdio.h"
#include "SuiteSparse_config.h"

namespace SuiteSparse_Mongoose
{

class Graph
{
public:
    /** CSparse3 Interoperability ********************************************/
    Int cs_n;                            /* # columns                        */
    Int cs_m;                            /* # rows                           */
    Int cs_nz;                           /* # triplet entries or -1          */
    Int cs_nzmax;                        /* max # nonzeros                   */

    /** Graph Data ***********************************************************/
    Int n;                               /* # vertices                       */
    Int nz;                              /* # edges                          */
    Int *p;                              /* Column pointers                  */
    Int *i;                              /* Row indices                      */
    Weight *x;                           /* Edge weight                      */
    Weight *w;                           /* Node weight                      */
    Weight X;                            /* Sum of edge weights              */
    Weight W;                            /* Sum of node weights              */

    Weight H;                            /* Heuristic max penalty to assess. */

    /** Partition Data *******************************************************/
    bool *partition;                     /* T/F denoting partition side      */
    Weight *vertexGains;                 /* Gains for each vertex            */
    Int *externalDegree;                 /* # edges lying across the cut     */
    Int *bhIndex;                        /* Index+1 of a vertex in the heap  */
    Int *bhHeap[2];                      /* Heap data structure organized by
                                            boundaryGains descending         */
    Int bhSize[2];                       /* Size of the boundary heap        */

    /** Cut Cost Metrics *****************************************************/
    Weight heuCost;                      /* cutCost + balance penalty        */
    Weight cutCost;                      /* Sum of edge weights in cut set   */
    Weight W0;                           /* Sum of partition 0 node weights  */
    Weight W1;                           /* Sum of partition 1 node weights  */
    Weight imbalance;                    /* Degree to which the partitioning
                                            is imbalanced, and this is
                                            computed as (0.5 - W0/W).        */

    /** Matching Data ********************************************************/
    Graph *parent;                       /* Link to the parent graph         */
    Int clevel;                          /* Coarsening level for this graph  */
    Int cn;                              /* # vertices in coarse graph       */
    Int *matching;                       /* Linked List of matched vertices  */
    Int *matchmap;                       /* Map from fine to coarse vertices */
    Int *invmatchmap;                    /* Map from coarse to fine vertices */
    Int *matchtype;                      /* Vertex's match classification
                                             0: Orphan
                                             1: Standard (random, hem, shem)
                                             2: Brotherly
                                             3: Community                    */

    /** Mark Data ************************************************************/
    Int *mark;                           /* O(n) mark array                  */
    Int markValue;                       /* Mark array can be cleared in O(k)
                                            by incrementing markValue.
                                            Implicitly, a mark value less than
                                            markValue is unmarked.           */

    /* Constructor & Destructor */
    void *operator new(size_t size, Graph* ptr){ return ptr; }
    Graph()
    {
        cs_n = cs_m = cs_nz = cs_nzmax = 0;
        n = nz = 0;
        p = NULL;
        i = NULL;
        x = NULL;
        w = NULL;
        X = 0.0;
        W = 0.0;
        H = 0.0;

        partition = NULL;
        vertexGains = NULL;
        externalDegree = NULL;
        bhIndex = NULL;
        bhHeap[0] = bhHeap[1] = NULL;
        bhSize[0] = bhSize[1] = 0;

        heuCost = 0.0;
        cutCost = 0.0;
        W0 = 0.0;
        W1 = 0.0;
        imbalance = 0.0;

        parent = NULL;
        clevel = 0;
        cn = 0;
        matching = NULL;
        matchmap = NULL;
        invmatchmap = NULL;
        matchtype = NULL;

        mark = NULL;
        markValue = 1;
    }

    static Graph *Create (
        const Int _n,
        const Int _nz
    )
    {
        Graph *ret = (Graph*) SuiteSparse_calloc(1, sizeof(Graph));
        if(!ret) return NULL;
        new (ret) Graph();

        int n =
        ret->n = 
        ret->cs_n = 
        ret->cs_m = _n;

        int nz =
        ret->nz = 
        ret->cs_nzmax = _nz;

        ret->cs_nz = -1; /* Compressed Column Format */
        ret->p = (Int*) SuiteSparse_malloc(n+1, sizeof(Int));
        ret->i = (Int*) SuiteSparse_malloc(nz, sizeof(Int));
        ret->x = (Weight*) SuiteSparse_malloc(nz, sizeof(Weight));
        ret->w = (Weight*) SuiteSparse_malloc(n, sizeof(Weight));
        ret->X = 0.0;
        ret->W = 0.0;
        ret->H = 0.0;
        if(!ret->p || !ret->i || !ret->x || !ret->w)
        {
             ret->~Graph();
             SuiteSparse_free(ret);
             return NULL;
        }

        ret->partition = (bool*) SuiteSparse_malloc(n, sizeof(bool));
        ret->vertexGains = (Weight*) SuiteSparse_malloc(n, sizeof(Weight));
        ret->externalDegree = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->bhIndex = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->bhHeap[0] = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->bhHeap[1] = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->bhSize[0] = ret->bhSize[1] = 0;
        if(!ret->partition || !ret->vertexGains || !ret->externalDegree
        || !ret->bhIndex || !ret->bhHeap[0] || !ret->bhHeap[1])
        {
             ret->~Graph();
             SuiteSparse_free(ret);
             return NULL;
        }

        ret->heuCost = 0.0;
        ret->cutCost = 0.0;
        ret->W0 = 0.0;
        ret->W1 = 0.0;
        ret->imbalance = 0.0;

        ret->parent = NULL;
        ret->clevel = 0;
        ret->cn = 0;
        ret->matching = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->matchmap = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->invmatchmap = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->matchtype = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->mark = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->markValue = 1;
        if(!ret->matching || !ret->matchmap || !ret->invmatchmap || !ret->mark
        || !ret->matchtype )
        {
             ret->~Graph();
             SuiteSparse_free(ret);
             return NULL;
        }

        return ret;
    }

    static Graph *Create (
        Graph *_parent
    )
    {
        Graph *ret = (Graph*) SuiteSparse_calloc(1, sizeof(Graph));
        if(!ret) return NULL;
        new (ret) Graph();

        ret->n = 
        ret->cs_n = 
        ret->cs_m = _parent->cn;
        int n = ret->n;

        ret->nz = 
        ret->cs_nzmax = _parent->nz;
        int nz = ret->nz;

        ret->cs_nz = -1; /* Compressed Column Format */
        ret->p = (Int*) SuiteSparse_malloc((n+1), sizeof(Int));
        ret->i = (Int*) SuiteSparse_malloc(nz, sizeof(Int));
        ret->x = (Weight*) SuiteSparse_malloc(nz, sizeof(Weight));
        ret->w = (Weight*) SuiteSparse_malloc(n, sizeof(Weight));
        ret->X = 0.0;
        ret->W = _parent->W;
        ret->H = 0.0;
        if(!ret->p || !ret->i || !ret->x || !ret->w)
        {
             ret->~Graph();
             SuiteSparse_free(ret);
             return NULL;
        }

        ret->partition = (bool*) SuiteSparse_malloc(n, sizeof(bool));
        ret->vertexGains = (Weight*) SuiteSparse_malloc(n, sizeof(Weight));
        ret->externalDegree = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->bhIndex = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->bhHeap[0] = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->bhHeap[1] = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->bhSize[0] = ret->bhSize[1] = 0;
        if(!ret->partition || !ret->vertexGains || !ret->externalDegree
        || !ret->bhIndex || !ret->bhHeap[0] || !ret->bhHeap[1])
        {
             ret->~Graph();
             SuiteSparse_free(ret);
             return NULL;
        }

        ret->heuCost = 0.0;
        ret->cutCost = 0.0;
        ret->W0 = 0.0;
        ret->W1 = 0.0;
        ret->imbalance = 0.0;

        ret->parent = _parent;
        ret->clevel = ret->parent->clevel+1;
        ret->cn = 0;
        ret->matching = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->matchmap = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->invmatchmap = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->matchtype = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->mark = (Int*) SuiteSparse_calloc(n, sizeof(Int));
        ret->markValue = 1;
        if(!ret->matching || !ret->matchmap || !ret->invmatchmap || !ret->mark
        || !ret->matchtype )
        {
             ret->~Graph();
             SuiteSparse_free(ret);
             return NULL;
        }

        return ret;
    }

    ~Graph()
    {
        this->n  = cs_n  = cs_m     = 0;
        this->nz = cs_nz = cs_nzmax = 0;

        p = (Int*) SuiteSparse_free(p);
        i = (Int*) SuiteSparse_free(i);
        x = (Weight*) SuiteSparse_free(x);
        w = (Weight*) SuiteSparse_free(w);

        partition = (bool*) SuiteSparse_free(partition);
        vertexGains = (Weight*) SuiteSparse_free(vertexGains);
        externalDegree = (Int*) SuiteSparse_free(externalDegree);
        bhIndex = (Int*) SuiteSparse_free(bhIndex);
        bhHeap[0] = (Int*) SuiteSparse_free(bhHeap[0]);
        bhHeap[1] = (Int*) SuiteSparse_free(bhHeap[1]);

        cutCost = 0.0;
        W0 = 0.0;
        W1 = 0.0;
        imbalance = 0.0;

        clevel = 0;
        cn = 0;
        matching = (Int*) SuiteSparse_free(matching);
        matchmap = (Int*) SuiteSparse_free(matchmap);
        invmatchmap = (Int*) SuiteSparse_free(invmatchmap);
        matchtype = (Int*) SuiteSparse_free(matchtype);

        mark = (Int*) SuiteSparse_free(mark);
        markValue = 0;
    }
};

}

/* Mongoose graph-related macros */
#ifndef MONGOOSE_MARKED
#define MONGOOSE_MARKED(a)   (mark[(a)] == markValue)
#endif

#ifndef MONGOOSE_MARK
#define MONGOOSE_MARK(a)     (mark[(a)] = markValue);
#endif

#endif
