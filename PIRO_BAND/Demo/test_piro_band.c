/* ========================================================================== */
/* === PIRO_BAND/Test/test_piro_band.c ====================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Simple test for band reduction of a small unsymmetric matrix.
 * */

#include "piro_band.h"
#include "piro_band_util.h"
#include "test_piro_band.h"

/*
 * Test the band reduction routine using a matrix stored in a file.
 * Input:
 *  sizefile - Size file name. The size file will have the dimensions of the
 *             the band matrix, the lower and upper bandwidths.
 *  pfile    - Matrix file name. The file with the numerical entries
 *             (real/complex) of the band matrix stored in packed band format.
 * */
static Int PIRO_BAND(test)()    /* returns 1 if OK, 0 if test failed */
{
    Int m, n, bl, bu ;
    Int ldab  ;
    COV_Entry *A , *work, *U, *V, *C, *Atemp ;
    Entry *D, *E, *diag ;
    COV_Entry *Aptr, *AllA ;
    COV_Entry *mtemp ;
    Int blk[4] ;
    Int dsize ;
    Int err ;

    m = 10 ;
    n = 10 ;
    bl = 4 ;
    bu = 4 ;
    printf("m = "ID", n= "ID", bl="ID", bu="ID"\n", m, n, bl, bu) ;

    ldab = bl+bu+1 ;

    /* Allocate memory two times the required memory for the matrix and
     * keep a copy as band reduction will destroy the input matrix. */
    AllA = (COV_Entry *) malloc( 2 * n * (ldab) * sizeof(COV_Entry)) ;
    if (!AllA)
    {
        printf("Out of memory") ;
        return 0 ;
    }

    Aptr = AllA ;
    A = Aptr ;
    Aptr += (n * ldab) ;
    Atemp = Aptr ;

    /* Get the numerical values from the file */
    if (!PIRO_BAND(get_random_matrix)(ldab, n, A, 0, bu))
    {
        printf("Error in generating matrix \n") ;
        free(AllA) ;
        return 0 ;
    }

    memcpy(Atemp, A, n * ldab * sizeof(COV_Entry)) ;

    /* Allocate space for the singular vectors */
    if (!PIRO_BAND(init_input)(m, n, &mtemp, &U, &V, &C))
    {
        printf("Error in initializing input \n") ;
        free(AllA) ;
        return 0 ;
    }

    /* Allocate memory bidiagonal (D, E). */
    diag = (Entry *) malloc( 2 * MIN(m, n) * sizeof(Entry)) ;
    if (!diag)
    {
        printf("Out of memory") ;
        free(AllA) ;
        free(mtemp) ;
    }
    D = diag ;
    E = diag + MIN(m, n) ;
    E[0] = 0.0 ;

    /* Get the suggested block size for the given matrix and allocate
     * memory for the workspace.
     * */
    PIRO_BAND_LONG_NAME(get_blocksize)(m, n, bl, bu, 1, blk) ;

    dsize = blk[0]*blk[1] > blk[2]*blk[3] ? blk[0]*blk[1] : blk[2]*blk[3] ;
    work = (COV_Entry *) malloc(2 * dsize * sizeof(COV_Entry)) ;
    if (!work)
    {
        free(AllA) ;
        free(mtemp) ;
        free(diag) ;
        printf("Out of memory") ;
        return 0 ;
    }

    /* Reduce the input matrix to the bidiagonal form */
    err = PIRO_BAND(reduce)(blk, m, n, m, bl, bu, (Entry *) A, ldab, D, E+1,
            (Entry *) U, m, (Entry *) V, n, (Entry *) C, m, (Entry *) work,
                0) ;
    if (err != 0)
    {
        printf("test Band reduction failed "ID"\n", err) ;
    }

    /* Check that the A - (U * S * V') is small enough */
    if (!PIRO_BAND(chk_output)(m, n, bu, bl, ldab, 0, D, E, U, V, C, Atemp, 0))
    {
        free(AllA) ;
        free(mtemp) ;
        free(diag) ;
        if (!work) free(work) ;
        printf("Chk_output failed\n") ;
        return 0 ;
    }

    /* Free Allocated space */
    free(AllA) ;
    free(mtemp) ;
    free(diag) ;
    if (!work) free(work) ;
    return 1 ;
}


/* Run the test for band reduction using a set of matrices. Runs a simple test.
 * See ../Tcov/ for the entire test coverage. */
Int PIRO_BAND(run_tests)( void )
{
    return PIRO_BAND(test)() ;
}
