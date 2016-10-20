#include "cs.h"
/* create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
 * where (I-beta*v*v')*x = s*e1 and e1 = [1 0 ... 0]'.
 * Note that this CXSparse version is different than CSparse.  See Higham,
 * Accuracy & Stability of Num Algorithms, 2nd ed, 2002, page 357. */
CS_ENTRY cs_house (CS_ENTRY *x, double *beta, CS_INT n)
{
    CS_ENTRY s;
    CS_INT i ;
    s = CS_ZERO();
    if (!x || !beta) return (CS_MINUS_ONE()) ;          /* check inputs */
    /* s = norm(x) */
    for (i = 0 ; i < n ; i++) s = CS_ADD(s, CS_MUL(x [i], CS_CONJ (x [i]))) ;
    s = CS_SQRT (s) ;
    if (CS_IS_ZERO(s))
    {
        (*beta) = 0 ;
        x [0] = CS_ONE() ;
    }
    else
    {
        /* s = sign(x[0]) * norm (x) ; */
        if (!CS_IS_ZERO(x [0]))
        {
            s = CS_MUL(s, CS_DIV(x [0], CS_MAKE_ENTRY(CS_ABS (x [0])))) ;
        }
        x [0] = CS_ADD(x [0], s) ;
        (*beta) = 1. / CS_REAL (CS_MUL(CS_CONJ (s), x [0])) ;
    }
    return (CS_MUL(CS_MINUS_ONE(), s)) ;
}
