#include "gp_mex.hpp"

namespace SuiteSparse_Mongoose
{

#define MEX_STRUCT_PUT(F)        addFieldWithValue(returner, #F, (Double) O->F);

mxArray *mex_put_options
(
    const Options *O
)
{
    mxArray *returner = mxCreateStructMatrix(1, 1, 0, NULL);

    MEX_STRUCT_PUT(randomSeed);
    MEX_STRUCT_PUT(coarsenLimit);
    MEX_STRUCT_PUT(matchingStrategy);
    MEX_STRUCT_PUT(doCommunityMatching);
//    MEX_STRUCT_PUT(davisBrotherlyThreshold);
    MEX_STRUCT_PUT(guessCutType);
    MEX_STRUCT_PUT(numDances);
    MEX_STRUCT_PUT(useFM);
    MEX_STRUCT_PUT(useQPGradProj);
//    MEX_STRUCT_PUT(gradprojTol);
//    MEX_STRUCT_PUT(gradprojIterationLimit);
//    MEX_STRUCT_PUT(targetSplit);
//    MEX_STRUCT_PUT(tolerance);

    return returner;
}

}