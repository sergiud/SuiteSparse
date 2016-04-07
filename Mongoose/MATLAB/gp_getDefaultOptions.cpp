
#include "gp_mex.hpp"

using namespace SuiteSparse_Mongoose;

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    cs Amatrix;
    int i, gtype, values;

    const char* usage = "Usage: O = gp_getDefaultOptions()";
    if(nargout != 1 || nargin != 0) mexErrMsgTxt(usage);

    Options *ret = Options::Create();
    if(ret == NULL)
        mexErrMsgTxt("Out of memory encountered while creating default Options");

    pargout[0] = mex_put_options(ret);
}
