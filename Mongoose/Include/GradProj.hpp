#ifndef GRADPROJ_HPP_
#define GRADPROJ_HPP_

#include "Mongoose_internal.hpp"
#include "QPDelta.hpp"

namespace SuiteSparse_Mongoose
{

void QPlinks
(
    Graph *G,
    Options *O,
    QPDelta *QP          /* pointer to QPDelta structure  */
);

Double QPgradproj
(
    Graph *G,
    Options *O,
    QPDelta *QP
);

void QPboundary
(
    Graph *G,
    Options *O,
    QPDelta *QP
);

#define QP_GRADPROJ_saveContextAndExit()                                      \
{                                                                             \
    QP->its = it;                                                             \
    QP->err = err;                                                            \
    QP->nf = nf;                                                              \
                                                                              \
    Double b = 0.;                                                            \
    if(ib != 0)                                                               \
    {                                                                         \
        b = (ib > 0 ? hi : lo);                                               \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        for (Int k = 0; k < n; k++) b += Ew[k] * x[k];                        \
    }                                                                         \
    QP->ib = ib;                                                              \
    QP->b = b;                                                                \
                                                                              \
    return err;                                                               \
}                                                                             \

}

#endif
