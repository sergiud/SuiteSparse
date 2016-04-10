#ifndef FIDDUCIAMATTHEYES_HPP_
#define FIDDUCIAMATTHEYES_HPP_

#include "Mongoose_internal.hpp"
#include "CutCost.hpp"

namespace SuiteSparse_Mongoose
{

/* Swap candidates have the following features: */
struct SwapCandidate
{
    Int vertex;
    bool partition;
    Weight nodeWeight;
    Weight gain;
    Weight heuCost;
    Int bhPosition;
    Weight imbalance;

    SwapCandidate(){
        vertex = 0;
        partition = false;
        nodeWeight = 0.0;
        gain = -INFINITY;
        heuCost = INFINITY;
        bhPosition = -1;
        imbalance = 0.0;
    }
};

void fmRefine
(
    Graph *G,
    Options *O
);

void fmSwap
(
    Graph *G,
    Options *O,
    Int vertex,
    Weight gain,
    bool oldPartition,
    Int *mark,
    Int markValue
);

void calculateGain
(
    Graph *G,
    Options *O,
    Int vertex,
    Weight *out_gain,
    Int *out_externalDegree
);

}

#endif /* FIDDUCIAMATTHEYES_HPP_ */
