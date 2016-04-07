#ifndef BOUNDARYHEAP_HPP_
#define BOUNDARYHEAP_HPP_

#include "Mongoose_internal.hpp"

/* Mongoose BoundaryHeap-related Macros */
#ifndef MONGOOSE_BOUNDARYHEAP_MACROS
#define MONGOOSE_BOUNDARYHEAP_MACROS

  #define MONGOOSE_HEAP_PARENT(a)    (((a)-1) / 2)
  #define MONGOOSE_LEFT_CHILD(a)     (2*(a) + 1)
  #define MONGOOSE_RIGHT_CHILD(a)    (2*(a) + 2)

  #define MONGOOSE_IN_BOUNDARY(v)     (bhIndex[(v)] > 0)
  #define MONGOOSE_PUT_BHINDEX(v,p)    bhIndex[(v)] = (p) + 1;
  #define MONGOOSE_GET_BHINDEX(v)     (bhIndex[(v)]-1)

#endif

namespace SuiteSparse_Mongoose
{

void bhLoad
(
    Graph *G,
    Options *O
);

void bhClear
(
    Graph *G
);

void bhInsert
(
    Graph *G,
    Int vertex
);

#if 0
Int bhDequeue
(
    Graph *G,
    Options *O,
    Int *bhIndex,
    Int *bhHeap,
    Int *inout_size,
    Weight *gains
);
#endif

void bhRemove
(
    Graph *G,
    Options *O,
    Int vertex,
    Weight gain,
    bool partition,
    Int bhPosition
);

void heapifyUp
(
    Int *bhIndex,
    Int *bhHeap,
    Weight *gains,
    Int vertex,
    Int position,
    Weight gain
);

void heapifyDown
(
    Int *bhIndex,
    Int *bhHeap,
    Int size,
    Weight *gains,
    Int vertex,
    Int position,
    Weight gain
);

}

#endif
