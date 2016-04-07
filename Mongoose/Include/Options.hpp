#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include "Mongoose_internal.hpp"

namespace SuiteSparse_Mongoose
{

struct Options
{
    Int randomSeed;

    /** Coarsening Options ***************************************************/
    Int coarsenLimit;
    MatchingStrategy matchingStrategy;
    bool doCommunityMatching;
    Weight davisBrotherlyThreshold;

    /** Guess Partitioning Options *******************************************/
    GuessCutType guessCutType;   /* The guess cut type to use */
    Int guessSearchDepth;        /* # of times to run pseudoperipheral node
                                    finder before settling on a guess        */

    /** Waterdance Options ***************************************************/
    Int numDances;               /* The number of interplays between FM and QP
                                    at any one coarsening level. */

    /**** Fidducia-Mattheyes Options *****************************************/
    bool useFM;                  /* Flag governing the use of FM             */
    Int fmSearchDepth;           /* The # of non-positive gain move to make  */
    Int fmConsiderCount;         /* The # of heap entries to consider        */
    Int fmMaxNumRefinements;     /* Max # of times to run FidduciaMattheyes  */

    /**** Quadratic Programming Options **************************************/
    bool useQPGradProj;          /* Flag governing the use of gradproj       */
    bool useQPBallOpt;           /* Flag governing the use of qp ball opt    */
    Weight gradprojTol;          /* Convergence tol for projected gradient   */
    Int gradprojIterationLimit;  /* Max # of iterations for gradproj         */

    /** Final Partition Target Metrics ***************************************/
    Weight targetSplit;          /* The desired split ratio (default 50/50)  */
    Weight tolerance;            /* The allowable tolerance (default 1%)     */

    /** Debug Options ********************************************************/
    bool doExpensiveChecks;

#if 0
    /** DOT Integration Options **********************************************/
    bool writeDot;
    const char *problemName;
    Int graphSID;
#endif

    /* Constructor & Destructor */
    static Options *Create()
    {
        Options *ret = (Options*) SuiteSparse_calloc(1, sizeof(Options));
        
        if(ret != NULL)
        {
            ret->randomSeed = 0;

            ret->coarsenLimit = 256;
            ret->matchingStrategy = HEMDavisPA;
            ret->doCommunityMatching = false;
            ret->davisBrotherlyThreshold = 2.0;

            ret->guessCutType = Pseudoperipheral_Fast;
            ret->guessSearchDepth = 10;

            ret->numDances = 1;

            ret->useFM = true;
            ret->fmSearchDepth = 50;
            ret->fmConsiderCount = 3;
            ret->fmMaxNumRefinements = 20;

            ret->useQPGradProj = true;
            ret->useQPBallOpt = true;
            ret->gradprojTol = 0.001;
            ret->gradprojIterationLimit = 50;

            ret->targetSplit = 0.5;
            ret->tolerance = 0.01;

            ret->doExpensiveChecks = false;
        }

        return ret;
    }
};

}

#endif
