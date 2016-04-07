/****************************************************************************
* metismex.c
* Public domain MATLAB CMEX-file to let you use METIS-4.0 from MATLAB.
* Usage:
* [part,edgecut] = metismex('PartGraphRecursive',A,nparts,wgtflag,options)
* [part,edgecut] = metismex('PartGraphKway',A,nparts,wgtflag,options)
* [perm,iperm] = metismex('EdgeND',A,options)
* [perm,iperm] = metismex('NodeND',A,options)
* sep = metismex('NodeBisect',A,wgtflag,options)
*
* Output arguments, along with the wgtflag and options input arguments,
* are optional. See the METIS manual for the meanings of the arguments.
*
* Note that error checking is not done: make sure A is structurally
* symmetric or it will crash.
*
* To compile, you need to have Metis 4, and do something like
*   mex -I<metis.h directory> -L<libmetis.a directory> -lmetis metismex.c
*
* Robert Bridson
*****************************************************************************/

#include "mex.h"
#include <strings.h>
#include "metis.h"

/*************************************************************************
* Given a graph, find a node separator bisecting it (roughly). The number
* of bisector nodes is returned in *nbnd, and the nodes themselves in the
* array bnds, which should already be allocated to have enough room.
**************************************************************************/
void METIS_NodeBisect(int nvtxs, idxtype *xadj, idxtype *adjncy,
                      idxtype *vwgt, idxtype *adjwgt, int wgtflag,
                      int *options, int *nbnd, idxtype *bnds, double unbalance)
{
  int i, j, tvwgt, tpwgts2[2];
  GraphType graph;
  CtrlType ctrl;
  idxtype *label, *bndind;

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType   = ONMETIS_CTYPE;
    ctrl.IType   = ONMETIS_ITYPE;
    ctrl.RType   = ONMETIS_RTYPE;
    ctrl.dbglvl  = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType   = options[OPTION_CTYPE];
    ctrl.IType   = options[OPTION_ITYPE];
    ctrl.RType   = options[OPTION_RTYPE];
    ctrl.dbglvl  = options[OPTION_DBGLVL];
  }
  ctrl.nseps = 5;    /* Take the best of 5 separators */
  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 50;

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  SetUpGraph(&graph, OP_ONMETIS, nvtxs, 1, xadj,adjncy,vwgt,adjwgt,wgtflag);

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph.vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  ctrl.maxvwgt = (1.5*tvwgt)/ctrl.CoarsenTo;

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, 2);

  MlevelNodeBisectionMultiple (&ctrl, &graph, tpwgts2, unbalance);

  IFSET(ctrl.dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph.nvtxs, graph.pwgts[0], graph.pwgts[1], graph.pwgts[2]));

  /* Now indicate the vertex separator */
  *nbnd = graph.nbnd;
  bndind = graph.bndind;
  label = graph.label;
  for (i = 0; i < *nbnd; ++i) 
    bnds[i] = label[bndind[i]];

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);
}

void convertMatrix (const mxArray *A, idxtype **xadj, idxtype **adjncy,
                    idxtype **vwgt, idxtype **adjwgt)
{
    int i, j, jbar, n, nnz, *jc, *ir;
    double *pr;

    /* Find MATLAB's matrix structure */
    n = mxGetN(A);
    jc = mxGetJc(A);
    nnz = jc[n];
    ir = mxGetIr(A);
    pr = mxGetPr(A);

    /* Allocate room for METIS's structure */
    *xadj = (idxtype*) mxCalloc (n+1, sizeof(idxtype));
    *adjncy = (idxtype*) mxCalloc (nnz, sizeof(idxtype));
    *vwgt = (idxtype*) mxCalloc (n, sizeof(idxtype));
    *adjwgt = (idxtype*) mxCalloc (nnz, sizeof(idxtype));

    /* Scan the matrix, not copying diagonal terms, and rounding doubles
     * to integer weights */
    (*xadj)[0] = 0;
    jbar = 0;
    for (i = 1; i <= n; i++) {
        for (j = jc[i-1]; j < jc[i]; j++) {
            if (ir[j] != i-1) {
                (*adjncy)[jbar] = ir[j];
                (*adjwgt)[jbar] = (int) pr[j];
                jbar++;
            } else {
                (*vwgt)[i-1] = (int) pr[j];
            }
        }
        (*xadj)[i] = jbar;
    }
}


#define FUNC_IN (prhs[0])
#define A_IN (prhs[1])
#define NPARTS_IN (prhs[2])
#define WGTFLAG_IN (prhs[3])
#define PARTOPTS_IN (prhs[4])
#define NDOPTS_IN (prhs[2])
#define NBWGTFLAG_IN (prhs[2])
#define NBOPTS_IN (prhs[3])

#define PART_OUT (plhs[0])
#define EDGECUT_OUT (plhs[1])
#define PERM_OUT (plhs[0])
#define IPERM_OUT (plhs[1])
#define SEP_OUT (plhs[0])

#define FUNCNAMELEN 25

/****************************************************************************
* mexFunction: gateway routine for MATLAB interface.
*****************************************************************************/
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, n, nparts, wgtflag, options[8] = {0}, numflag = 0, edgecut, sepsize;
    idxtype *xadj, *adjncy, *vwgt, *adjwgt, *part, *perm, *iperm, *sep;
    char funcname[FUNCNAMELEN];
    double *optarray, *partpr, *permpr, *ipermpr, *seppr;

    /* First do some general argument checking */
    if (nrhs < 2 || nrhs > 5 || nlhs > 2) {
        mexErrMsgTxt ("Wrong # of arguments");
    }
    if (!mxIsChar(FUNC_IN)) {
        mexErrMsgTxt ("First parameter must be a string");
    }
    n = mxGetN(A_IN);
    if (!mxIsSparse(A_IN) || n!=mxGetM(A_IN)) {
        mexErrMsgTxt ("Second parameter must be a symmetric sparse matrix");
    }

    /* Copy the matrix over, getting rid of diagonal, and converting to
     * integer weights */
    convertMatrix (A_IN, &xadj, &adjncy, &vwgt, &adjwgt);

    /* Now figure out which function we have to do */
    mxGetString (FUNC_IN, funcname, FUNCNAMELEN);

    if (strcasecmp(funcname,"PartGraphRecursive")==0) {

        /* Figure out values for nparts, wgtflag and options */
        if (nrhs < 3) {
            mexErrMsgTxt ("Third parameter needed: nparts");
        }
        nparts = (int) mxGetScalar (NPARTS_IN);
        if (nrhs >= 4) {
            wgtflag = (int) mxGetScalar (WGTFLAG_IN);
        } else {
            wgtflag = 0;
        }
        if (nrhs >= 5) {
            optarray = mxGetPr (PARTOPTS_IN);
            for (i = 1; i < 4; ++i) {
                options[i] = (int) optarray[i-1];
            }
        }

        /* Allocate memory for result of call */
        part = (idxtype*) mxCalloc (n, sizeof(idxtype));

        /* Do the call */
        METIS_PartGraphRecursive (&n, xadj, adjncy, vwgt, adjwgt, &wgtflag,
                                  &numflag, &nparts, options, &edgecut, part);

        /* Figure out output values */
        if (nlhs >= 1) {
            PART_OUT = mxCreateDoubleMatrix (1, n, mxREAL);
            partpr = mxGetPr (PART_OUT);
            for (i = 0; i < n; i++) {
                partpr[i] = (double) part[i];
            }

            if (nlhs >= 2) {
                EDGECUT_OUT = mxCreateDoubleMatrix (1, 1, mxREAL);
                mxGetPr(EDGECUT_OUT)[0] = (double) edgecut;
            }
        }

    } else if (strcasecmp(funcname,"PartGraphKway")==0) {

        /* Figure out values for nparts, wgtflag and options */
        if (nrhs < 3) {
            mexErrMsgTxt ("Third parameter needed: nparts");
        }
        nparts = (int) mxGetScalar (NPARTS_IN);
        if (nrhs >= 4) {
            wgtflag = (int) mxGetScalar (WGTFLAG_IN);
        } else {
            wgtflag = 0;
        }
        if (nrhs >= 5) {
            optarray = mxGetPr (PARTOPTS_IN);
            for (i = 1; i < 4; ++i) {
                options[i] = (int) optarray[i-1];
            }
        }

        /* Allocate memory for result of call */
        part = (idxtype*) mxCalloc (n, sizeof(idxtype));

        /* Do the call */
        METIS_PartGraphKway (&n, xadj, adjncy, vwgt, adjwgt, &wgtflag,
                             &numflag, &nparts, options, &edgecut, part);

        /* Figure out output values */
        if (nlhs >= 1) {
            PART_OUT = mxCreateDoubleMatrix (1, n, mxREAL);
            partpr = mxGetPr (PART_OUT);
            for (i = 0; i < n; i++) {
                partpr[i] = (double) part[i];
            }

            if (nlhs >= 2) {
                EDGECUT_OUT = mxCreateDoubleMatrix (1, 1, mxREAL);
                mxGetPr(EDGECUT_OUT)[0] = (double) edgecut;
            }
        }

    } else if (strcasecmp(funcname,"EdgeND")==0) {

        /* Figure out values for options */
        if (nrhs >= 3) {
            optarray = mxGetPr (NDOPTS_IN);
            for (i = 1; i < 4; ++i) {
                options[i] = (int) optarray[i-1];
            }
        }

        /* Allocate memory for result of call */
        perm = (idxtype*) mxCalloc (n, sizeof(idxtype));
        iperm = (idxtype*) mxCalloc (n, sizeof(idxtype));

        /* Do the call */
        METIS_EdgeND (&n, xadj, adjncy, &numflag, options, perm, iperm);

        /* Figure out output values */
        if (nlhs >= 1) {
            PERM_OUT = mxCreateDoubleMatrix (1, n, mxREAL);
            permpr = mxGetPr (PERM_OUT);
            for (i = 0; i < n; i++) {
                permpr[i] = perm[i]+1.0;
            }

            if (nlhs >= 2) {
                IPERM_OUT = mxCreateDoubleMatrix (1, n, mxREAL);
                ipermpr = mxGetPr (IPERM_OUT);
                for (i = 0; i < n; i++) {
                    ipermpr[i] = iperm[i]+1.0;
                }
            }
        }

    } else if (strcasecmp(funcname,"NodeND")==0) {

        /* Figure out values for options */
        if (nrhs >= 3) {
            optarray = mxGetPr (NDOPTS_IN);
            for (i = 1; i < 4; ++i) {
                options[i] = (int) optarray[i-1];
            }
            for (i = 5; i < 8; ++i) {
                options[i] = (int) optarray[i-2];
            }
        }

        /* Allocate memory for result of call */
        perm = (idxtype*) mxCalloc (n, sizeof(idxtype));
        iperm = (idxtype*) mxCalloc (n, sizeof(idxtype));

        /* Do the call */
        METIS_NodeND (&n, xadj, adjncy, &numflag, options, perm, iperm);

        /* Figure out output values */
        if (nlhs >= 1) {
            PERM_OUT = mxCreateDoubleMatrix (1, n, mxREAL);
            permpr = mxGetPr (PERM_OUT);
            for (i = 0; i < n; i++) {
                permpr[i] = perm[i]+1.0;
            }

            if (nlhs >= 2) {
                IPERM_OUT = mxCreateDoubleMatrix (1, n, mxREAL);
                ipermpr = mxGetPr (IPERM_OUT);
                for (i = 0; i < n; i++) {
                    ipermpr[i] = iperm[i]+1.0;
                }
            }
        }

    } else if (strcasecmp(funcname,"NodeBisect")==0) {

        if (nrhs >= 3) {
            wgtflag = (int) mxGetScalar (NBWGTFLAG_IN);
        } else {
            wgtflag = 0;
        }
        if (nrhs >= 4) {
            optarray = mxGetPr (NBOPTS_IN);
            for (i = 1; i < 4; ++i) {
                options[i] = (int) optarray[i-1];
            }
        }

        /* Allocate memory for result of call */
        sep = (idxtype*) mxCalloc (n, sizeof(idxtype));

        /* Do the call */
        METIS_NodeBisect (n, xadj, adjncy, vwgt, adjwgt, wgtflag,
                          options, &sepsize, sep, 1.5);

        /* Figure out output values */
        if (nlhs >= 1) {
            SEP_OUT = mxCreateDoubleMatrix (1, sepsize, mxREAL);
            seppr = mxGetPr (PART_OUT);
            for (i = 0; i < sepsize; i++) {
                seppr[i] = (double) (sep[i]+1);
            }
        }

    } else {
        mexErrMsgTxt ("Unknown metismex function");
    }
}

