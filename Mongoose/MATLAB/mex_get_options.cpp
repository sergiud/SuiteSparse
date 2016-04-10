#include "gp_mex.hpp"

namespace SuiteSparse_Mongoose
{

#define MEX_STRUCT_READINT(F)    returner->F = (Int) readField(matOptions, #F);
#define MEX_STRUCT_READDOUBLE(F) returner->F = readField(matOptions, #F);
#define MEX_STRUCT_READBOOL(F)   returner->F = (readField(matOptions, #F) != 0.0);

Options *mex_get_options
(
    const mxArray *matOptions
)
{
    Options *returner = Options::Create();
    if(!returner)
        return NULL;
    if(matOptions == NULL)
        return returner;

    MEX_STRUCT_READINT(randomSeed);
    MEX_STRUCT_READINT(coarsenLimit);
    returner->matchingStrategy = (MatchingStrategy) readField(matOptions, "matchingStrategy");
    MEX_STRUCT_READBOOL(doCommunityMatching);
//    MEX_STRUCT_READDOUBLE(davisBrotherlyThreshold);
    returner->guessCutType = (GuessCutType) readField(matOptions, "guessCutType");
    MEX_STRUCT_READINT(numDances);
    MEX_STRUCT_READBOOL(useFM);
    MEX_STRUCT_READBOOL(useQPGradProj);
//    MEX_STRUCT_READDOUBLE(gradprojTol);
//    MEX_STRUCT_READINT(gradprojIterationLimit);
//    MEX_STRUCT_READDOUBLE(targetSplit);
//    MEX_STRUCT_READDOUBLE(tolerance);

    return returner;
}

}