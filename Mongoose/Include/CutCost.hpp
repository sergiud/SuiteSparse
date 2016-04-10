#ifndef CUTCOST_HPP_
#define CUTCOST_HPP_

namespace SuiteSparse_Mongoose
{

/* A partitioning has the following metrics: */
struct CutCost
{
    Weight heuCost;     /* Sum of cutCost and balance penalty     */
    Weight cutCost;     /* Sum of edge weights in the cut set.    */
    Weight W[2];        /* Sum of node weights in each partition. */
    Weight imbalance;   /* targetSplit - (W[0] / W)               */
};

}

#endif /* CUTCOST_HPP_ */
