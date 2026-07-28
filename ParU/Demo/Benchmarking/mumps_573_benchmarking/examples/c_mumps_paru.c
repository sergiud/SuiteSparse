//superlu_mt : ourbranch : pdlinsolx.c
//suitesparse:dev2 branch: paru/demo/paru_benchmark.cpp : line 453

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "omp.h"
#include "mmio.h"
#include "dmumps_c.h"

#define NTRIALS 3
#define THREAD_CONFIGS 6
//#define THREAD_CONFIGS 7

#define USE_COMM_WORLD -987654
#define CHECK_ERROR( id , phase )\
{\
    if ( id.infog[0] < 0 )\
    {\
        printf( " Error in the %s phase\tINFOG(1)= %d, INFOG(2)= %d\n", phase , id.infog[0] , id.infog[1] ) ;\
        exit( 1 ) ;\
    }\
}

// Comparison function for qsort (for ascending order)
int compare(const void *a, const void *b)
{
    return ( *(double *) a - *(double *) b ) ;
}

// Function to find the median of an array
double median( double arr[], int n )
{
    // Sort the array in ascending order
    qsort( arr , n , sizeof(double) , compare ) ;

    if ( n % 2 == 0 ) // If the number of elements is even
    {
        // Median is the average of the two middle elements
        return ( arr [ n / 2 - 1 ] + arr[ n / 2 ] ) / 2.0 ;
    }
    else // If the number of elements is odd
    {
        // Median is the middle element
        return arr[ n / 2 ] ;
    }
}

// This routine is copied / modified from https://math.nist.gov/MatrixMarket/mmio/c/example_read.c
int read_matrix
(
 char *file_name ,
 DMUMPS_STRUC_C *id
)
{
    /*
     * Read Matrix Market file
     */

    MUMPS_INT m , n ;
    MUMPS_INT8 nnz ;
    MUMPS_INT *irn ;
    MUMPS_INT *jcn ;
    double *a ;

    FILE *f ;
    int ret_code ;
    MM_typecode matcode ;
    MUMPS_INT8 i ;

    if ( ( f = fopen( file_name , "r" ) ) == NULL )
    {
       exit( 1 ) ;
    }

    if ( mm_read_banner( f , &matcode ) != 0 )
    {
        printf("Could not process Matrix Market banner.\n") ;
        exit( 1 ) ;
    }

    if ( mm_is_complex( matcode ) )
    {
        printf( "Sorry, only real matrices are supported." ) ;
        printf( "Matrix Market type: [%s]\n", mm_typecode_to_str( matcode ) ) ;
        exit( 1 ) ;
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nnz)) !=0)
        exit(1);
    
    if ( m != n )
    {
        printf( "Matrix must be square ( m = %d , n = %d )\n" , m , n ) ;
        exit( 1 ) ;
    }

    /* reseve memory for matrix */

    irn = (MUMPS_INT *) malloc(nnz * sizeof(MUMPS_INT));
    jcn = (MUMPS_INT *) malloc(nnz * sizeof(MUMPS_INT));
    a = (double *) malloc(nnz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nnz; i++)
    {
        fscanf(f, "%d %d %lg\n", &irn[i], &jcn[i], &a[i]) ;
        irn[i] ;
        jcn[i] ;
    }

    if (f !=stdin) fclose(f);

    id->n = n ;
    id->nnz = nnz ;
    id->irn = irn ;
    id->jcn = jcn ;
    id->a = a ;
}

int generate_rhs
(
 DMUMPS_STRUC_C *id
)
{
    id->rhs = ( double * ) malloc ( id->n * sizeof( double ) ) ;
    if ( ! id->rhs )
    {
        printf ( "Not enough memory for RHS allocation.\n" ) ;
        exit ( 1 ) ;
    }

    // initialize the right-hand-side
    for ( MUMPS_INT i = 0 ; i < id->n ; i++ )
    {
        id->rhs[ i ] = i + 1 ;
    }

    return 0 ;
}

#if defined(MAIN_COMP)
/*
 * Some Fortran compilers (COMPAQ fort) define "main" in
 * their runtime library while a Fortran program translates
 * to MAIN_ or MAIN__ which is then called from "main".
 * We defined argc/argv arbitrarily in that case.
 */
int MAIN__();
int MAIN_()
    {
        return MAIN__();
    }

int MAIN__()
{
    int argc=1;
    char * name = "c_mumps_paru";
    char ** argv ;
#else
int main(int argc, char ** argv)
{
#endif

    int nth[THREAD_CONFIGS] = { 1 , 2 , 4 , 8 , 16 , 24 } ;
//    int nth[THREAD_CONFIGS] = { 1 , 2 , 5 , 10 , 20 , 40 , 80 } ;
    double trial_times[NTRIALS] ;
    double ana_time ;
    double fac_times[THREAD_CONFIGS] ;
    double sol_times[THREAD_CONFIGS] ;
    double tic , toc ;

    DMUMPS_STRUC_C id ;
    int error ;

    // Check input

    if ( argc < 2 )
    {
        fprintf( stderr , "Usage: %s [martix-market-filename]\n", argv[0] ) ;
        exit( 1 ) ;
    }

    /* MPI initialize */

#if defined(MAIN_COMP)
    argv = &name;
#endif
    error = MPI_Init(&argc, &argv);

    int rank ;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
    if ( rank > 0 )
    {
        printf( "This program must be run sequentially, with a single MPI process.\n" ) ;
        exit( 1 ) ;
    }

    /* MUMPS initialize */

    id.comm_fortran = USE_COMM_WORLD;
    id.par = 1 ;
    id.sym = 0 ;
    id.job = -1 ;
    dmumps_c( &id ) ;
    CHECK_ERROR( id , "Initialize" ) ;

    /* Problem initialize */

    error = read_matrix( argv[1] , &id ) ;
    error = generate_rhs( &id ) ;

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    // Centralized assembled matrix
    id.ICNTL(5) = 0 ; id.ICNTL(18) = 0 ;
    /* No outputs */
//  id.ICNTL(1) = -1 ; id.ICNTL(2) = -1 ; id.ICNTL(3) = -1 ; id.ICNTL(4) = 0 ;
    /* Verbose messages */
    id.ICNTL(4) = 3 ;
    // compute statistics
    id.ICNTL(11) = 2 ;

    int ordering_code_in_MUMPS[ 2 ] = { 0 , 5 } ;
    // ordering [ 0 ] = 0 => AMD
    // ordering [ 1 ] = 5 => METIS
    for ( int ordering_code_in_PARU = 0 ; ordering_code_in_PARU < 2 ; ordering_code_in_PARU++ )
    {
        for ( int nthid = 0 ; nthid < THREAD_CONFIGS ; nthid++ )
        {
            printf ("TESTING MUMPS: ORDERING (ParU code): %d, THREADS: %d\n",
                ordering_code_in_PARU, nth [ nthid ] ) ;

            omp_set_num_threads( nth [ nthid ] ) ;

            /* MUMPS analyse */

            for ( int trial = 0 ; trial < NTRIALS ; trial++ )
            {
                id.ICNTL(7) = ordering_code_in_MUMPS [ ordering_code_in_PARU ] ;
                // Activate multi-threading
                id.ICNTL(48) = 1 ;
                id.ICNTL(16) = nth [ nthid ] ;
                id.job = 1 ;
                tic = omp_get_wtime( ) ;
                dmumps_c( &id ) ;
                toc = omp_get_wtime( ) ;
                trial_times[ trial ] = toc - tic ;
                CHECK_ERROR( id , "Analysis" ) ;
            }
            ana_time = median( trial_times , NTRIALS ) ;

            /* MUMPS factorization */

            for ( int trial = 0 ; trial < NTRIALS ; trial++ )
            {
                id.job = 2 ;
                tic = omp_get_wtime( ) ;
                dmumps_c( &id ) ;
                toc = omp_get_wtime( ) ;
                trial_times[ trial ] = toc - tic ;
                CHECK_ERROR( id , "Factorization" ) ;
            }
            fac_times[ nthid ] = median( trial_times , NTRIALS ) ;

            /* MUMPS solve*/

            for ( int trial = 0 ; trial < NTRIALS ; trial++ )
            {

                // initialize the right-hand-side
                for ( MUMPS_INT i = 0 ; i < id.n ; i++ )
                {
                    id.rhs[ i ] = i + 1 ;
                }

                id.job = 3 ;
                tic = omp_get_wtime( ) ;
                dmumps_c( &id ) ;
                toc = omp_get_wtime( ) ;
                trial_times[ trial ] = toc - tic ;
                CHECK_ERROR( id , "Solve" ) ;
                printf ("scaled resid: %g\n", id.rinfog [5]) ;
            }
            sol_times[ nthid ] = median( trial_times , NTRIALS ) ;
        }

        /* Print benchmark timings */

        printf ("TABLE, MUMPS, %s, %d, sym_time:, %12.6e,",
            (argv[1]== NULL) ? " " : argv[1],
            ordering_code_in_PARU, ana_time) ;
        printf (" num_times:, ") ;
        for (int kk = 0 ; kk < THREAD_CONFIGS ; kk++)
        {
            if (fac_times [kk] < 0) break ;
            printf (" %12.6e, ", fac_times [kk]) ;
        }
        printf (" sol_times:, ") ;
        for (int kk = 0 ; kk < THREAD_CONFIGS ; kk++)
        {
            if (sol_times [kk] < 0) break ;
            printf (" %12.6e, ", sol_times [kk]) ;
        }
        printf ("\n") ;
    }

    /* Problem terminate */

    free( id.irn ) ;
    free( id.jcn ) ;
    free( id.a   ) ;
    free( id.rhs   ) ;

    /* MUMPS terminate */

    id.job = -2 ;
    dmumps_c( &id ) ;
    CHECK_ERROR( id , "Terminate" ) ;

    /* MPI finalize */

    error = MPI_Finalize() ;

    return 0 ;
}
