
#include "mex.h"
#include "Mongoose.hpp"

namespace SuiteSparse_Mongoose
{

void shcpDataToMAT
(
    mxArray* matStruct,
    const char* field,
    mxClassID classID,
    void* data,
    size_t size
);

void addFieldWithValue
(
    mxArray* matStruct,     /* The mxArray assumed to be a matlab structure. */
    const char* fieldname,  /* The name of the field to create.              */
    const double value      /* The double value to assign to the new field.  */
);

double readField
(
    const mxArray* matStruct,
    const char* fieldname
);

Options *mex_get_options(const mxArray *Omatlab = NULL);
mxArray *mex_put_options(const Options *O);

Graph *mex_get_graph
(
    const mxArray *Gmatlab,        /* The sparse matrix            */
    const mxArray *Amatlab = NULL  /* The real-valued node weights */
);

Int *gp_mex_get_int
(
    Int n,
    const mxArray *Imatlab,
    Int *imax,
    Int lo
);

mxArray *gp_mex_put_int(Int *p, Int n, Int offset, Int do_free);
mxArray *gp_mex_put_logical(bool *p, Int n);
double *gp_mex_get_double (Int n, const mxArray *X) ;
double *gp_mex_put_double (Int n, const double *b, mxArray **X) ;

}