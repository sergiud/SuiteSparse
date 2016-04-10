#ifndef GUESSCUT_HPP_
#define GUESSCUT_HPP_

#include "Mongoose_internal.hpp"
#include "CutCost.hpp"
#include "BoundaryHeap.hpp"

namespace SuiteSparse_Mongoose
{

bool guessCut(Graph *G, Options *O);

void pseudoperipheralGuess(Graph *G, Options *O);
void findAllPseudoperipheralNodes(Graph *G, Options *O, Int *list, Int *listsize);

Int diagBFS
(
    Graph *G,
    Options *O,
    Int *stack,
    Int *mark,
    Int markStart,
    Int *start
);

void partBFS
(
    Graph *G,
    Options *O,
    Int start
);

}

#endif
