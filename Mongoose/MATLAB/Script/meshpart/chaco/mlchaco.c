/* MLCHACO : Mex-file version of Hendrickson/Leland's graph partitioner.
 *
 * [map,t]=mlchaco(A, vwgts, ewgtsP, xy, assignment, architecture, ...
 *   ndims_tot, mesh_dims, goal, global_method, local_method, rqi_flag, ...
 *   vmax, ndims, eigtol, seed);
 *
 * A:          Sparse adjacency matrix of graph (symmetric, zero diagonal).
 * vwgts:      If non-null, vertex weights.
 * ewgtsP:     Nonzero means use elements of A as edge weights.
 * xy:         Each row of xy is the coordinates of a vertex.
 * map:        Output partition, as a dense map vector.
 *
 * Other arguments are as described in Section 6.10 of the Chaco manual.  
 * A must be sparse; all other args must be dense matrices or scalars.
 *
 * This routine is meant to be called from the Matlab routine "chaco.m",
 * but a user can call it directly if he or she really wants to.  The
 * arguments are almost the same as the arguments to Chaco's "interface"
 * routine, which this routine calls.
 *
 * This is part of a Matlab interface to the software described in
 * B. Hendrickson and R. Leland, "The Chaco User's Guide (Version 2.0)",
 * Sandia National Laboratories report SAND94-2692, October 1994.
 * Interface written by John Gilbert, Xerox PARC, February 1995, and
 * Copyright (c) 1994-1995 by Xerox Corporation.  All rights reserved.
 * See copyright.m for complete copyright and licensing notice.
 *
 * Modified by Tim Davis, July 1999.  Added 2nd output argument,
 * which is the cpu time taken by chaco itself.
 */

#include <stdio.h>
#include "mex.h"
#include <sys/times.h>

/* Aliases for input and output arguments */
#define A_in               prhs[0]
#define vwgts_in           prhs[1]
#define ewgtsP_in          prhs[2]
#define xy_in              prhs[3]
#define assignment_in      prhs[4]
#define architecture_in    prhs[5]
#define ndims_tot_in       prhs[6]
#define mesh_dims_in       prhs[7]
#define goal_in            prhs[8]
#define global_method_in   prhs[9]
#define local_method_in    prhs[10]
#define rqi_flag_in        prhs[11]
#define vmax_in            prhs[12]
#define ndims_in           prhs[13]
#define eigtol_in          prhs[14]
#define seed_in            prhs[15]
#define map_out            plhs[0]


void mexFunction(    
    int         nlhs,           /* number of expected outputs */
    Matrix      *plhs[],        /* matrix pointer array returning outputs */
    int         nrhs,           /* number of inputs */
    Matrix      *prhs[]         /* matrix pointer array for inputs */
    )
{
    /* Arguments to Chaco "interface" routine: */

    int       nvtxs;          /* number of vertices in full graph */
    int      *start;          /* start of edge list for each vertex */
    int      *adjacency;      /* edge list data */
    int      *vwgts;          /* weights for all vertices */
    float    *ewgts;          /* weights for all edges */
    float    *x, *y, *z;      /* coordinates for inertial method */
    char     *outassignname;  /* name of assignment output file */
    char     *outfilename;    /* output file name */
    short    *assignment;     /* set number of each vtx (length n) */
    int       architecture;   /* 0 => hypercube, d => d-dimensional mesh */
    int       ndims_tot;      /* total number of cube dimensions to divide */
    int       mesh_dims[3];   /* dimensions of mesh of processors */
    double   *goal;           /* desired set sizes for each set */
    int       global_method;  /* global partitioning algorithm */
    int       local_method;   /* local partitioning algorithm */
    int       rqi_flag;       /* should I use RQI/Symmlq eigensolver? */
    int       vmax;           /* if so, how many vertices to coarsen down to? */
    int       ndims;          /* number of eigenvectors (2^d sets) */
    double    eigtol;         /* tolerance on eigenvectors */
    long      seed;           /* for random graph mutations */
 
    int    i, m, nedges, geodims, failure;
    double *p;

    double *Tout ;			/* time: added by Tim Davis */
    struct tms t1 ;
    struct tms t2 ;

    /* Check for correct number of arguments passed from Matlab. */
    if (nrhs != 16) {
        mexErrMsgTxt("MLCHACO requires 16 input arguments.");
    } else if (nlhs > 2) {
        mexErrMsgTxt("MLCHACO requires one or two output arguments.");
    }

    /* Slog through the arguments, converting from Matlab matrices
       to the various types Chaco wants.                            */

    nvtxs = mxGetN(A_in);
    start = mxGetJc(A_in);
    nedges = start[nvtxs];
    adjacency = mxGetIr(A_in);
    if (mxGetN(vwgts_in)) {
        vwgts = (int *) mxCalloc(nvtxs, sizeof(int));
        p = mxGetPr(vwgts_in);
        for (i = 0; i < nvtxs; vwgts[i++] = (int) *p++);
    } else vwgts = NULL;
    if (mxGetScalar(ewgtsP_in)) {
        ewgts = (float *) mxCalloc(nedges, sizeof(float));
        p = mxGetPr(A_in);
        for (i = 0; i < nedges; ewgts[i++] = (float) *p++);
    } else ewgts = NULL;
    geodims = mxGetN(xy_in);
    if (geodims >= 1) {
        x = (float *) mxCalloc(nvtxs, sizeof(float));
        p = mxGetPr(xy_in);
        for (i = 0; i < nvtxs; x[i++] = (float) *p++);
    } else x = NULL;
    if (geodims >= 2) {
        y = (float *) mxCalloc(nvtxs, sizeof(float));
        for (i = 0; i < nvtxs; y[i++] = (float) *p++);
    } else y = NULL;
    if (geodims >= 3) {
        z = (float *) mxCalloc(nvtxs, sizeof(float));
        for (i = 0; i < nvtxs; z[i++] = (float) *p++);
    } else z = NULL;
    outassignname = NULL;
    outfilename = NULL;
    assignment = (short *) mxCalloc(nvtxs, sizeof(short));
    if (mxGetN(assignment_in)) {
        p = mxGetPr(assignment_in);
        for (i = 0; i < nvtxs; assignment[i++] = (short) *p++);
    }
    architecture = (int) mxGetScalar(architecture_in);
    ndims_tot = (int) mxGetScalar(ndims_tot_in);
    p = mxGetPr(mesh_dims_in);
    m = mxGetM(mesh_dims_in);
    if (m < mxGetN(mesh_dims_in))
        m = mxGetN(mesh_dims_in);
    for (i = 0; i < m; mesh_dims[i++] = (int) *p++);
    if (mxGetN(goal_in)) {
        goal = mxGetPr(goal_in);
    } else goal = NULL;
    global_method = (int) mxGetScalar(global_method_in);
    local_method = (int) mxGetScalar(local_method_in);
    rqi_flag = (int) mxGetScalar(rqi_flag_in);
    vmax = (int) mxGetScalar(vmax_in);
    ndims = (int) mxGetScalar(ndims_in);
    eigtol = (double) mxGetScalar(eigtol_in);
    seed = (long) mxGetScalar(seed_in);

    /* Finally, call Chaco */

    (void) times (&t1) ;

    failure = interface(nvtxs, start, adjacency, vwgts, ewgts, x, y, z,
        outassignname, outfilename, assignment, architecture, ndims_tot, 
        mesh_dims, goal, global_method, local_method, rqi_flag, vmax, 
        ndims, eigtol, seed);

    (void) times (&t2) ;

    if (nlhs > 1)
    {
	plhs [1] = mxCreateFull (1, 1, REAL) ;
	Tout = mxGetPr (plhs [1]) ;
	Tout [0] =
	((double) (t2.tms_utime + t2.tms_stime - t1.tms_utime - t1.tms_stime)) /
	((double) CLK_TCK) ;
    }

    /* Copy the mapping vector into a Matlab matrix. */

    if (!failure) {
        map_out = mxCreateFull(1,nvtxs,REAL);
        p = mxGetPr(map_out);
        for (i = 0; i < nvtxs; *p++ = (double) assignment[i++]);
    }

    /* Free what we allocated */
   
    if (vwgts != NULL) mxFree((char *) vwgts);
    if (ewgts != NULL) mxFree((char *) ewgts);
    if (x != NULL) mxFree((char *) x);
    if (y != NULL) mxFree((char *) y);
    if (z != NULL) mxFree((char *) z);
    if (assignment != NULL) mxFree((char *) assignment);

    if (failure)
        mexErrMsgTxt("The call to C Chaco returned failure.");

    return;
}
