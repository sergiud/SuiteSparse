#ifndef INTEROP_HPP_
#define INTEROP_HPP_

#include "string.h"
#include "cs.hpp"
#include "Mongoose_internal.hpp"

namespace SuiteSparse_Mongoose
{

/* Take a matrix market format file, sanitize it, and use the sanitized file
 * to load a CSparse3 compressed column sparse matrix and returns it. */
Graph *readGraphFromMM(const char *mmFilename);

/* load a 1-indexed triplet matrix from a flat text file */
cs *cs_load2 (FILE *f);

/* Configure an on-stack CSparse matrix from an existing Mongoose Graph. */
cs *GraphToCSparse3(Graph *G, bool copy = false);
Graph *CSparse3ToGraph(cs *G, bool resetEW = false, bool resetNW = false);

#if 0
/* DOT is the language used in sfdp, fdp, dot, etc. */
enum DOTRenderSetting{ GraphOnly, CutSet, Matching, Heaps };
void writeDot
(
    Graph *G,
    Options *O,
    const char *comment,
    DOTRenderSetting in_setting = GraphOnly
);

void heapWriteDot
(
    Graph *G,
    Options *O,
    Int *heap,
    Int size,
    Weight *gains,
    Int vHighlight = -1
);
#endif

}

#endif
