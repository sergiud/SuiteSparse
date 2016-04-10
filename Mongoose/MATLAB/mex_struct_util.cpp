#include "gp_mex.hpp"

namespace SuiteSparse_Mongoose
{

/*****************************************************************************
Function:
    Shallow Copy Data to MATLAB

Purpose:
    Shallow-copies data into a new matlab structure field.

Assumes:
    mxREAL data only
    M x N is really just 1 x size
    the matlab structure is a single cell.
 *****************************************************************************/
void shcpDataToMAT
(
    mxArray* matStruct,
    const char* field,
    mxClassID classID,
    void* data,
    size_t size
)
{
    mxArray *newField;

    /* If we weren't given a matlab structure, do nothing. */
    if(!mxIsStruct(matStruct)) return;

    newField = mxCreateNumericMatrix(0, 0, classID, mxREAL);

    /* Brain Transplant */
    mxFree(mxGetData(newField));
    mxSetData(newField, data);
    mxSetM(newField, 1);
    mxSetN(newField, size);
    mexMakeMemoryPersistent(data);

    /* Add the new field to the given matlab structure. */
    mxAddField(matStruct, field);
    mxSetField(matStruct, 0, field, newField);
}

/*****************************************************************************
Function:
    Add Field with Value

Purpose:
    Adds a field to the given mxArray (assumed to be a 1x1 matlab structure)
    and automatically sets it to the given value.
 *****************************************************************************/
void addFieldWithValue
(
    mxArray* matStruct,     /* The mxArray assumed to be a matlab structure. */
    const char* fieldname,  /* The name of the field to create.              */
    const double value      /* The double value to assign to the new field.  */
)
{
    if(!mxIsStruct(matStruct)) return;

    mxAddField(matStruct, fieldname);
    mxSetField(matStruct, 0, fieldname, mxCreateDoubleScalar(value));
}

/*****************************************************************************
Function:
    Reads Field

Purpose:
    Reads the scalar value out of the specified structure.field.
 *****************************************************************************/
double readField
(
    const mxArray* matStruct,
    const char* fieldname
)
{
    double returner = 0.0;

    mxArray* field = mxGetField(matStruct, 0, fieldname);
    if(field != NULL) returner = mxGetScalar(field);
    return returner;
}

}