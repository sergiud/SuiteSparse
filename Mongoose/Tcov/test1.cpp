
#include "Mongoose.hpp"
#include "Conditioning.hpp"
#include "SuiteSparse_config.h"

using namespace SuiteSparse_Mongoose;

void RunAllTests (
    const char *inputFile, 
    Options *O, 
    int allowedMallocs
);

void RunTest (
    const char *inputFile, 
    Options *O,
    int allowedMallocs
);

/* Custom memory management functions allow for memory testing. */
int AllowedMallocs;

void *myMalloc(size_t size)
{
    if(AllowedMallocs <= 0) return NULL;
    AllowedMallocs--;
    return malloc(size);
}

void *myCalloc(size_t count, size_t size)
{
    if(AllowedMallocs <= 0) return NULL;
    AllowedMallocs--;
    return calloc(count, size);
}

void *myRealloc(void *ptr, size_t newSize)
{
    if(AllowedMallocs <= 0) return NULL;
    AllowedMallocs--;
    return realloc(ptr, newSize);
}

void myFree(void *ptr)
{
    if(ptr != NULL) free(ptr);
}

int main(int argn, const char **argv)
{
    SuiteSparse_start();

    /* Get the input file from the console input. */
    const char *inputFile = NULL;
    if(argn == 2) { inputFile = argv[1]; }
    else { printf("usage: Mongoose MM-input-file\n"); return 0; }    

    Options *O = Options::Create();
    if(!O)
    {
         SuiteSparse_finish();
         return 1; // Return an error if we failed.
    }

    O->doExpensiveChecks = true;

    /* Override SuiteSparse memory management with custom testers. */
    SuiteSparse_config.malloc_func = myMalloc;
    SuiteSparse_config.calloc_func = myCalloc;
    SuiteSparse_config.realloc_func = myRealloc;
    SuiteSparse_config.free_func = myFree;

    for(int i=0; i<200; i++)
    {
        RunAllTests(inputFile, O, i);
    }

    O->~Options();
    SuiteSparse_free(O);

    /* Return success */
    SuiteSparse_finish();
    return 0;
}

void RunAllTests (
    const char *inputFile,
    Options *O,
    int allowedMallocs
)
{
    printf("Running all tests with %d AllowedMallocs\n", allowedMallocs);

   /* 
    * (1 TEST)
    * TRY A BOGUS LOAD
    */
    RunTest("bogus", O, allowedMallocs);

    /* 
     * (12 TESTS)
     * TRY A VARIETY OF MATCHING STRATEGIES AND GUESS CUT STRATEGIES
     */
    MatchingStrategy matchingStrategies[4] = {Random, HEM, HEMPA, HEMDavisPA};
    GuessCutType guessCutStrategies[3] = {Pseudoperipheral_All, Pseudoperipheral_Fast, QP_GradProj};

    for(int c=0; c<2; c++)
    {
        O->doCommunityMatching = c;

        for(int i=0; i<4; i++)
        {
            O->matchingStrategy = matchingStrategies[i];

            for(int j=0; j<3; j++)
            {
                O->guessCutType = guessCutStrategies[j];

                printf("+ Testing %d %d %d\n", c, i, j);
                RunTest(inputFile, O, allowedMallocs);
                printf("- Testing %d %d %d\n", c, i, j);
            }
        }
    }

   /* 
    * (1 TEST)
    * DO A DEEP COPY OF THE GRAPH AND RESET EDGEWEIGHTS
    */
    Graph *DCG = conditionGraph(readGraphFromMM(inputFile), O, true);
    if(DCG)
    {
        cs *csDCG = GraphToCSparse3(DCG, true);
        csDCG = cs_spfree(csDCG);
        DCG->~Graph();
        SuiteSparse_free(DCG);
    }
}

void RunTest (
    const char *inputFile, 
    Options *O,
    int allowedMallocs
)
{
    /* Set the # of mallocs that we're allowed. */
    AllowedMallocs = allowedMallocs;

    /* Read and condition the matrix from the MM file. */
    Graph *inG = readGraphFromMM(inputFile);
    if(!inG) return;

    Graph *U = conditionGraph(inG, O);
    if(!U) return;

    ComputeEdgeSeparator(U, O);

    printf("U->cutCost = %f\n", U->cutCost);
    printf("U->imbalance = %f\n", U->imbalance);
    printf("Partitioning Complete\n");

    U->~Graph();
    SuiteSparse_free(U);
}
