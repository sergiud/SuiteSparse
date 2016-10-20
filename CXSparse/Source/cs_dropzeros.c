#include "cs.h"
static CS_INT cs_nonzero (CS_INT i, CS_INT j, CS_ENTRY aij, void *other)
{
    return (!CS_IS_ZERO(aij)) ;
}
CS_INT cs_dropzeros (cs *A)
{
    return (cs_fkeep (A, &cs_nonzero, NULL)) ;  /* keep all nonzero entries */
} 
