#ifndef HAGER_HPP_
#define HAGER_HPP_

#include "Mongoose_internal.hpp"
#include "QPDelta.hpp"
#include "GradProj.hpp"

namespace SuiteSparse_Mongoose
{

void qpGradProj
(
    Graph *G,
    Options *O,
    bool isInitial = false
);

void qpBallOpt
(
    Graph *G,
    Options *O
);

}

#endif
