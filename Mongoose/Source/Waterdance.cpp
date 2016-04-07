
#include "Waterdance.hpp"

namespace SuiteSparse_Mongoose
{

void waterdance(Graph *G, Options *O)
{
    Int numDances = O->numDances;
    for(Int i=0; i<numDances; i++)
    {
        fmRefine(G, O);
        qpGradProj(G, O);
    }
}

}