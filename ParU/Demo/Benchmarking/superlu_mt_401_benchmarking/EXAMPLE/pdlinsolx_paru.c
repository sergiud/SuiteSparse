/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
#include "slu_mt_ddefs.h"

// Revised by Tim Davis, from EXAMPLE/pdlinsolx.c:
#include "omp.h"
extern void dreadrb (int_t *, int_t *, int_t *, double **, int_t **, int_t **);
#ifndef _OPENMP
#error "OpenMP is required"
#endif

#include <unistd.h>
#include <math.h>

static int compar (const void *p1, const void *p2)
{
    double x1 = *((double *) p1) ;
    double x2 = *((double *) p2) ;
    return (x1 < x2 ? -1 : ((x1 > x2) ? 1 : 0)) ;
}

void
parse_command_line(int argc, char *argv[], int_t *nprocs, int_t *lwork, 
		   int_t *w, int_t *relax, double *u, fact_t *fact, 
		   trans_t *trans, yes_no_t *refact, equed_t *equed) ;

int
main(int argc, char *argv[])
{
    SuperMatrix A, L, U;
    SuperMatrix B, X;
    NCformat    *Astore;
    SCPformat   *Lstore;
    NCPformat   *Ustore;
    int_t         nprocs;
    fact_t      fact;
    trans_t     trans;
    yes_no_t    refact, usepr;
    equed_t     equed;
    double      *a;
    int_t         *asub, *xa;
    int_t         *perm_c; /* column permutation vector */
    int_t         *perm_r; /* row permutations from partial pivoting */
    void        *work;
    superlumt_options_t superlumt_options;
    int_t         info, lwork, nrhs, ldx, panel_size, relax;
    int_t         m, n, nnz, permc_spec, i;
    double      *rhsb, *rhsx, *xact;
    double      *R, *C;
    double      *ferr, *berr;
    double      u, drop_tol, rpg, rcond;
    superlu_memusage_t superlu_memusage;

    /* Default parameters to control factorization. */
    nprocs = 1;
    fact  = EQUILIBRATE;
    trans = NOTRANS;
    equed = NOEQUIL;
    refact= NO;
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    u     = 1.0;
    usepr = NO;
    drop_tol = 0.0;
    lwork = 0;
    nrhs  = 1;

    /* Command line options to modify default behavior. */
    parse_command_line(argc, argv, &nprocs, &lwork, &panel_size, &relax, 
		       &u, &fact, &trans, &refact, &equed);

    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
	if ( !work ) {
	    SUPERLU_ABORT("DLINSOLX: cannot allocate work[]");
	}
    }

#if ( PRNTlevel==1 )
    cpp_defs();
    printf("int_t %d bytes\n", sizeof(int_t));
#endif

#define RB
#if defined( DEN )
    m = n;
    nnz = n * n;
    dband(n, n, nnz, &a, &asub, &xa);
#elif defined( BAND )
    m = n;
    nnz = (2*b+1) * n;
    dband(n, b, nnz, &a, &asub, &xa);
#elif defined( BD )
    nb = 5;
    bs = 200;
    m = n = bs * nb;
    nnz = bs * bs * nb;
    dblockdiag(nb, bs, nnz, &a, &asub, &xa);
#elif defined( HB )
    dreadhb(&m, &n, &nnz, &a, &asub, &xa);
#elif defined( RB )
    dreadrb(&m, &n, &nnz, &a, &asub, &xa);
#else    
    dreadmt(&m, &n, &nnz, &a, &asub, &xa);
#endif

    // make a copy of the input matrix values
    double *Asave = doubleMalloc (nnz) ;
    memcpy (Asave, a, nnz * sizeof (double)) ;

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = A.Store;
    printf("Dimension " IFMT "x" IFMT "; # nonzeros " IFMT "\n", A.nrow, A.ncol, Astore->nnz);

    double *Aval = Astore->nzval ;

/*
    {
        printf ("\nA matrix as read in:\n") ;
        for (int j = 0 ; j < A.ncol ; j++)
        {
            for (int p = Astore->colptr [j] ; p < Astore->colptr [j+1] ; p++)
            {
                int_t i = Astore->rowind [p] ;
                printf ("%d %d %32.16e\n", (int) i, (int) j, Aval [p]) ;
            }
        }
    }
*/

    if (!(rhsb = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsb[].");
    if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);
    /*
    for (int k = 0 ; k < n ; k++)
    {
        printf ("xact [%d] = %g\n", k, xact [k]) ;
    }
    */

    // make a copy of the input right-hand-side
    DNformat *Bstore = B.Store ;
    double *Bmat = Bstore->nzval ;
    double *Bsave = doubleMalloc (m) ;
    memcpy (Bsave, Bmat, m * sizeof (double)) ;

    if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
    if (!(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double)))) 
        SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
        SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) ) 
        SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");

double resid = 0 ;
int_t nprocs_max = omp_get_max_threads ( ) ;
printf ("Max # of threads: %d\n", nprocs_max) ;
for (permc_spec = 3 ; permc_spec <= 3 ; permc_spec++)
{

double superlu_tanalyze_all [20] ;
double superlu_tfactor_all  [20] ;
double superlu_tsolve_all   [20] ;
int kthreads = 0 ;
for (int kk =  0 ; kk < 20 ; kk++)
{
    superlu_tanalyze_all [kk] = -1 ;
    superlu_tfactor_all  [kk] = -1 ;
    superlu_tsolve_all   [kk] = -1 ;
}

for (int_t nprocs = nprocs_max ; nprocs > 0 ; nprocs = (nprocs == 24) ? 16 : (nprocs/2))
{
printf ("\n============ SuperLU # of threads: %d Ordering: %d\n", nprocs, permc_spec) ;
#define NTRIALS 3
int middle = NTRIALS / 2 ;

double superlu_tanalyze [NTRIALS] ;
double superlu_tfactor  [NTRIALS] ;
double superlu_tsolve   [NTRIALS] ;

for (int trial = 0 ; trial < NTRIALS ; trial++)
{
printf ("\n=== Trial: %d\n", trial) ;

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
//  permc_spec = 1;
    double t1 = omp_get_wtime ( ) ;
    get_perm_c(permc_spec, &A, perm_c);

    superlumt_options.nprocs = nprocs;
    superlumt_options.fact = fact;
    superlumt_options.trans = trans;
    superlumt_options.refact = refact;
    superlumt_options.panel_size = panel_size;
    superlumt_options.relax = relax;
    superlumt_options.usepr = usepr;
    superlumt_options.drop_tol = drop_tol;
    superlumt_options.diag_pivot_thresh = u;
    superlumt_options.SymmetricMode = NO;
    superlumt_options.PrintStat = NO;
    superlumt_options.perm_c = perm_c;
    superlumt_options.perm_r = perm_r;
    superlumt_options.work = work;
    superlumt_options.lwork = lwork;
    if ( !(superlumt_options.etree = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for etree[].");
    if ( !(superlumt_options.colcnt_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    if ( !(superlumt_options.part_super_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    t1 = omp_get_wtime ( ) - t1 ;
    printf ("Order time =   %16.8f\n", t1) ;
    rcond = 0 ;
    // double t2 = omp_get_wtime ( ) ;

    DNformat *Xstore = X.Store ;
    double *Xmat = Xstore->nzval ;
    double *Resid = xact ;

    /* 
     * Solve the system and compute the condition number
     * and error bounds using pdgssvx.
     */
    pdgssvx(nprocs, &superlumt_options, &A, perm_c, perm_r,
	    &equed, R, C, &L, &U, &B, &X, &rpg, &rcond,
	    ferr, berr, &superlu_memusage, &info);

    superlu_tanalyze [trial] = t1 ;
    superlu_tfactor  [trial] = paru_benchmark_timings [0] ;
    superlu_tsolve   [trial] = paru_benchmark_timings [1] ;

    printf ("equed: %d, info: %d\n", equed, info) ;

    // restore the values of A and B; they are overwritten by pdgssvx
    // with the equilbriated matrix R*A*C and equilbriated B !!
    memcpy (Aval, Asave, nnz * sizeof (double)) ;
    memcpy (Bmat, Bsave, m * sizeof (double)) ;

    if ( info == 0 || info == n+1 ) {

	printf("Recip. pivot growth = %e\n", rpg);
	printf("Recip. condition number = %e\n", rcond);
	printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
	for (i = 0; i < nrhs; ++i) {
	    printf(IFMT "%16e%16e\n", i+1, ferr[i], berr[i]);
	}
	       
        Lstore = (SCPformat *) L.Store;
        Ustore = (NCPformat *) U.Store;
	printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
	       superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
	       superlu_memusage.expansions);
	     
	fflush(stdout);

        // compute the 1-norm of the residual, Resid = norm1 (B-A*X),
        // using xact as workspace

        // Resid = B
        // memcpy (Resid, Bmat, A.nrow * sizeof (double)) ;
        for (int k = 0 ; k < A.nrow ; k++)
        {
            Resid [k] = Bmat [k] ;
        }
        // Resid = Resid - A*X
        sp_dgemv ("N", (double) -1.0, &A, Xmat, (int_t) 1, (double) 1.0,
            Resid, (int_t) 1) ;
        resid = 0 ;
        double xnorm = 0 ;
        for (int k = 0 ; k < A.nrow ; k++)
        {
            resid += fabs (Resid [k]) ;
            xnorm += fabs (Xmat [k]) ;
            // if (k < 10)
            //  printf ("X (%d+1) = %g ;\n", k, Xmat [k]) ;
        }
//      for (int k = 0 ; k < A.nrow ; k++)
//          printf ("B (%d+1) = %g ;\n", k, Bmat [k]) ;
        printf ("absolute resid: %g (norm1 (b-Ax))\n", resid) ;
        printf ("norm1(x):       %g\n", xnorm);
        double anorm = 0 ;
        // printf ("\nA matrix:\n") ;
        for (int j = 0 ; j < A.ncol ; j++)
        {
            double colsum = 0 ;
            for (int p = Astore->colptr [j] ; p < Astore->colptr [j+1] ; p++)
            {
                int_t i = Astore->rowind [p] ;
                colsum += fabs (Aval [p]) ;
                // if (p < 100)
                //     printf ("%d %d %32.16e\n", (int) i, (int) j, Aval [p]) ;
            }
            if (colsum > anorm) { anorm = colsum ; }
        }
        printf ("norm1(A):       %g\n", anorm) ;
        resid = resid / (anorm * xnorm) ;
        printf ("Relative resid: %g (norm1 (b-Ax) / norm1(A)*norm1(x))\n\n", resid) ;

    } else if ( info > 0 && lwork == -1 ) {
        printf("pdgssvx(): info " IFMT "\n", info);
        printf("** Estimated memory: " IFMT " bytes\n", info - n);
    }

    SUPERLU_FREE (superlumt_options.etree);
    SUPERLU_FREE (superlumt_options.colcnt_h);
    SUPERLU_FREE (superlumt_options.part_super_h);
    if ( lwork == 0 ) {
        Destroy_SuperNode_SCP(&L);
        Destroy_CompCol_NCP(&U);
    }
}

    if ( info == 0 || info == n+1 ) {
        qsort (superlu_tanalyze, NTRIALS, sizeof (double), compar) ;
        qsort (superlu_tfactor , NTRIALS, sizeof (double), compar) ;
        qsort (superlu_tsolve  , NTRIALS, sizeof (double), compar) ;
        printf ("Aznaveh: th: %2d ord %d ", nprocs, permc_spec) ;

        switch (permc_spec)
        {
            case 0: printf ("NAT    ") ; break ; // natural ordering 
            case 1: printf ("MMDATA ") ; break ; // minimum degree ordering on structure of A'*A
            case 2: printf ("MMD    ") ; break ; //  minimum degree ordering on structure of A'+A
            case 3: printf ("COLAMD ") ; break ; // approximate minimum degree for unsymmetric matrices
        }

        printf ("analyze: %16.8f fac: %16.8f sol: %16.8f relresid %8.3e\n",
            superlu_tanalyze [middle],
            superlu_tfactor  [middle],
            superlu_tsolve   [middle],
            resid) ;

        superlu_tanalyze_all [kthreads] = superlu_tanalyze [middle] ;
        superlu_tfactor_all  [kthreads] = superlu_tfactor  [middle] ;
        superlu_tsolve_all   [kthreads] = superlu_tsolve   [middle] ;
    }
    kthreads++ ;

}

    printf ("TABLE, SuperLU_MT, threads:, %2d, ordering:, %2d, ", nprocs_max, permc_spec) ;

    printf (" analyze_time:, ") ;
    for (int kk = 0 ; kk < 20 ; kk++)
    {
        if (superlu_tanalyze_all [kk] < 0) break ;
        printf (" %12.8e, ", superlu_tanalyze_all [kk]) ;
    }

    printf (" factor_time:, ") ;
    for (int kk = 0 ; kk < 20 ; kk++)
    {
        if (superlu_tfactor_all [kk] < 0) break ;
        printf (" %12.8e, ", superlu_tfactor_all [kk]) ;
    }

    printf (" solve_time:, ") ;
    for (int kk = 0 ; kk < 20 ; kk++)
    {
        if (superlu_tsolve_all [kk] < 0) break ;
        printf (" %12.8e, ", superlu_tsolve_all [kk]) ;
    }
    printf ("\n") ;


}


    SUPERLU_FREE (Asave);
    SUPERLU_FREE (Bsave);
    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }
}

/*  
 * Parse command line options.
 */
void
parse_command_line(int argc, char *argv[], int_t *nprocs, int_t *lwork, 
		   int_t *w, int_t *relax, double *u, fact_t *fact, 
		   trans_t *trans, yes_no_t *refact, equed_t *equed)
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "hp:l:w:x:u:f:t:r:e:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-p <int> - number of processes\n");
	    printf("\t-l <int> - length of work[*] array\n");
	    printf("\t-w <int> - panel size\n");
	    printf("\t-x <int> - maximum size of relaxed supernodes\n");
	    printf("\t-u <int> - pivoting threshold\n");
	    printf("\t-f <FACTORED/DOFACT/EQUILIBRATE> - factor control\n");
	    printf("\t-t <NOTRANS/TRANS/CONJ> - transpose or not\n");
	    printf("\t-r <NO/YES> - refactor or not\n");
	    printf("\t-e <NOEQUIL/ROW/COL/BOTH> - equilibrate or not\n");
	    exit(1);
	    break;
	  case 'p': *nprocs = atoi(optarg);
	            break;
	  case 'l': *lwork = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'x': *relax = atoi(optarg); 
	            break;
	  case 'u': *u = atof(optarg); 
	            break;
	  case 'f': *fact = (fact_t) atoi(optarg);
	            break;
	  case 't': *trans = (trans_t) atoi(optarg);
	            break;
	  case 'r': *refact = (yes_no_t) atoi(optarg);
	            break;
	  case 'e': *equed = (equed_t) atoi(optarg);
	            break;
	  default: fprintf(stderr, "Invalid command line option.\n");
		   break;
  	}
    }
}
