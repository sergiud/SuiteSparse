/* ========================================================================== */
/* === QPnapsack2 =========================================================== */
/* ==========================================================================
    Solve a pair of projection problems:

        min ||x_0 - y_0|| subject to 0 <= x_0 <= 1,
                                     a_0'x_0 = b_0, lo_0 <= b_0 <= hi_0

        min ||x_1 - y_1|| subject to 0 <= x_1 <= 1,
                                     a_1'x_1 = b_1, lo_1 <= b_1 <= hi_1

    where x_0 stands for the first m components of an n component vector
    while x_1 stands for the last n-m components.  It is assumed that
    the column vector a is strictly positive since, in our application,
    the vector a represents the node weights, which are >= 1. If a is NULL,
    then it is assumed that a is identically 1.
   ========================================================================== */

#include "Napsack.hpp"

namespace SuiteSparse_Mongoose
{

void QPnapsack2
(
    Double *x,         /* holds y on input, and the solution x on output */
    Double *lambda,    /* initial guess (input) final value (output) multiplier
                          lambda [0], lambda [1] */
    Int *ix,           /* ix_i = +1, -1, or 0 on input, x_i =1, 0, or 0 < x_i < 1 */
    Int n,             /* # cols ? */
    Int m,             /* # rows ? */
    Double *a,         /* edge weights ? */
    Int *lo,           /* Problem lo limits in 0, 1 */
    Int *hi,           /* Problem hi limits in 0, 1 */
    Double *wx1,       /* Double workspace #1 */
    Int *wi1,          /* Integer workspace #1 */
    Int *wi2           /* Integer workspace #2 */
)
{
    Double *am;
    Int *ixm;

    lambda[0] = QPnapsack(x, m, lo[0], hi[0], a, lambda[0], ix, wx1, wi1, wi2);
    if (m >= n) return;

    am = (a ? a+m : NULL);
    ixm = (ix ? ix+m : NULL);

    lambda[1] = QPnapsack(x+m, n-m, lo[1], hi[1], am, lambda[1], ixm, wx1, wi1, wi2);
}

}
