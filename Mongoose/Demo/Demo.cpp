
#include "Mongoose.hpp"
#include "Conditioning.hpp"

using namespace SuiteSparse_Mongoose;

int main(int argn, const char **argv)
{
    /* Get the input file from the console input. */
    const char *inputFile = NULL;
    if(argn == 2) { inputFile = argv[1]; }
    else { printf("usage: Mongoose MM-input-file\n"); return 0; }

    Options *O = Options::Create();
    if(!O) return 1; // Return an error if we failed.

    O->doExpensiveChecks = true;
    O->matchingStrategy = HEMDavisPA;
    O->guessCutType = QP_BallOpt;

    /* Read & condition the matrix market input. */
    Graph *U = conditionGraph(readGraphFromMM(inputFile), O);

    ComputeEdgeSeparator(U, O);

printf("U->cutCost = %f\n", U->cutCost);
printf("U->imbalance = %f\n", U->imbalance);
printf("Partitioning Complete\n");

    U->~Graph();
    free(U);

    O->~Options();
    free(O);

    /* Return success */
    return 0;
}
