/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Modified by John Gilbert 25 Aug 1995 to change "exit" call to        *
 * mex-file error message.                                              */

#include	<stdio.h>
#include	<string.h>
#include "mex.h"

/* Wrapper for exit() - print message and exit with status code. Exit code
   of 0 indicates normal termination. Exit code of 1 indicates early 
   termination following detection of some problem. Call with bail(NULL,status) 
   to suppress message. */ 
void      bail(msg, status)
char     *msg;
int       status;
{
    extern FILE *Output_File;		/* Output file or NULL */

    if (msg != NULL && (int) strlen(msg) > 0) {
	if (Output_File != NULL) {
            fprintf(Output_File, "%s\n", msg);
	}
    }
    mexErrMsgTxt(msg);
}
