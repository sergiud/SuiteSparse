#include "cs.h"
/* apply the ith Householder vector to x */
CS_INT cs_happly (const cs *V, CS_INT i, double beta, CS_ENTRY *x)
{
    CS_INT p, *Vp, *Vi ;
    CS_ENTRY *Vx, tau;
    tau = CS_ZERO();
    if (!CS_CSC (V) || !x) return (0) ;     /* check inputs */
    Vp = V->p ; Vi = V->i ; Vx = V->x ;
    for (p = Vp [i] ; p < Vp [i+1] ; p++)   /* tau = v'*x */
    {
        tau = CS_ADD(tau, CS_MUL(CS_CONJ (Vx [p]), x [Vi [p]])) ;
    }
    tau = CS_MUL(tau, CS_MAKE_ENTRY(beta)) ;                           /* tau = beta*(v'*x) */
    for (p = Vp [i] ; p < Vp [i+1] ; p++)   /* x = x - v*tau */
    {
        x [Vi [p]] = CS_SUB(x[Vi[p]], CS_MUL(Vx [p], tau)) ;
    }
    return (1) ;
}
