#ifndef QPDelta_HPP_
#define QPDelta_HPP_

#include <new>
#include "SuiteSparse_config.h"
using namespace SuiteSparse_Mongoose;

namespace SuiteSparse_Mongoose
{

class QPDelta
{
public:
    Double *x;                 /* current estimate of solution            */
    Int nf;                    /* number of i such that 0 < x_i < 1       */
    Int *ix;                   /* ix_i = +1,-1, or 0 if x_i = 1,0, or 0 < x_i < 1 */
    Int *LinkUp;               /* linked list for free indices            */
    Int *LinkDn;               /* linked list, LinkDn [LinkUp [i]] = i    */
    Double *g;                 /* gradient at current x                   */
    Double *D;                 /* max value along the column.             */

    Int *wi[2];
    Double *wx[3];

    Int its;
    Double err;
    Int ib;
    Double b;

    static QPDelta *Create(Int n)
    {
        QPDelta *ret = (QPDelta*) SuiteSparse_calloc(1, sizeof(QPDelta));
        if(!ret) return NULL;

        ret->x = (Double*) SuiteSparse_malloc(n, sizeof(Double));
        ret->ix = (Int*) SuiteSparse_malloc(n, sizeof(Int));
        ret->LinkUp = (Int*) SuiteSparse_malloc(n+1, sizeof(Int));
        ret->LinkDn = (Int*) SuiteSparse_malloc(n+1, sizeof(Int));
        ret->g = (Double*) SuiteSparse_malloc(n, sizeof(Double));
        ret->D = (Double*) SuiteSparse_malloc(n, sizeof(Double));
        ret->wi[0] = (Int*) SuiteSparse_malloc(n+1, sizeof(Int));
        ret->wi[1] = (Int*) SuiteSparse_malloc(n+1, sizeof(Int));
        ret->wx[0] = (Double*) SuiteSparse_malloc(n, sizeof(Double));
        ret->wx[1] = (Double*) SuiteSparse_malloc(n, sizeof(Double));
        ret->wx[2] = (Double*) SuiteSparse_malloc(n, sizeof(Double));

        if(!ret->x || !ret->ix || !ret->LinkUp || !ret->LinkDn 
        || !ret->g || !ret->D || !ret->wi[0] || !ret->wi[1] 
        || !ret->wx[0] || !ret->wx[1] || !ret->wx[2])        
        {
            ret->~QPDelta();
            ret = (QPDelta*) SuiteSparse_free(ret);
        }

        return ret;   
    }

    ~QPDelta()
    {
        x = (Double*) SuiteSparse_free(x);
        ix = (Int*) SuiteSparse_free(ix);
        LinkUp = (Int*) SuiteSparse_free(LinkUp);
        LinkDn = (Int*) SuiteSparse_free(LinkDn);
        g = (Double*) SuiteSparse_free(g);
        D = (Double*) SuiteSparse_free(D);

        for(Int i=0; i<2; i++)
            wi[i] = (Int*) SuiteSparse_free(wi[i]);

        for(Int i=0; i<3; i++)
            wx[i] = (Double*) SuiteSparse_free(wx[i]);
    }
};

}

#endif
