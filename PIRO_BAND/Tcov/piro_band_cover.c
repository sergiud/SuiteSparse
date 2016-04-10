/* ========================================================================== */
/* === PIRO_BAND/Tcov/test_piro_band.c ====================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */


#include "piro_band.h"
#include "piro_band_util.h"
#include "piro_band_lapack.h"
#include "piro_band_lapack_internal.h"
#include "piro_band_cover.h"
#include "test_piro_band.h"

/* ========================================================================== */
/* === piro_band_test ======================================================= */
/* ========================================================================== */

/*
 * Test the band reduction routine using a matrix stored in a file.
 * Input:
 *  m, n, bl, bu:  matrix is m-by-n with lower and upper bandwidths bl and bu.
 *  sym      - 0/1. A flag that is one if the input matrix is symmetric and
 *             zero if it is unsymmetric.
 *  testmem  - 0/1. A flag to test the MALLOC in piro_band_reduce.
 *  tblk     - block size to override the recommended block size. If NULL
 *             use recommended block size.
 * */
static Int PIRO_BAND(test)  /* return 1 if OK, 0 if failure */
(
    Int m,
    Int n,
    Int bl,
    Int bu,
    Int sym,            /* create a symmetric test matrix */
    Int testmem,        /* test memory allocation */
    Int *tblk,
    int *ntests
)
{
    Int ldab  ;
    COV_Entry *A , *work, *U, *V, *C, *T ;
    Entry *D, *E, *diag ;
    COV_Entry *Aptr, *AllA ;
    COV_Entry *mtemp ;
    Int blk[4] ;
    Int *iblk ;
    Int dsize ;
    Int err ;

    printf("m = "ID", n= "ID", bl="ID", bu="ID"\n", m, n, bl, bu) ;
    (*ntests) ++ ;

    /* Always use the upper triangular part for the symmetric case. */
    if (sym)
    {
        bl = 0 ;
    }

    ldab = bl+bu+1 ;

    /* Allocate memory two times the required memory for the matrix and
     * keep a copy as band reduction will destroy the input matrix. */
    AllA = (COV_Entry *) malloc( 2 * n * (ldab) * sizeof(COV_Entry)) ;
    if (!AllA)
    {
        printf("Out of memory") ;
        return (0) ;
    }

    Aptr = AllA ;
    A = Aptr ;
    Aptr += (n * ldab) ;
    T = Aptr ;

    /* Get the numerical values from the file */
    if (!PIRO_BAND(get_random_matrix) (ldab, n, A, sym, bu))
    {
        printf("Error in generating matrix \n") ;
        free(AllA) ;
        return (0) ;
    }

    memcpy(T, A, n * ldab * sizeof(COV_Entry)) ;

    /* Allocate space for the singular vectors */
    if (!PIRO_BAND(init_input) (m, n, &mtemp, &U, &V, &C))
    {
        printf("Error in initializing input \n") ;
        free(AllA) ;
        return (0) ;
    }

    /* Allocate memory for the workspace, bidiagonal (D, E). */
    diag = (Entry *) malloc( 2 * MIN(m, n) * sizeof(Entry)) ;
    if (!diag)
    {
        printf("Out of memory") ;
        free(AllA) ;
        free(mtemp) ;
        return (0) ;
    }
    D = diag ;
    E = diag + MIN(m, n) ;
    E[0] = 0.0 ;

    if (!testmem && tblk == NULL)
    {
        /* Get the suggested block size for the given matrix  */
        /* Always use wantuv = 0, as wantuv = 1 will be covered by smaller
         * matrices anyway. */
        PIRO_BAND_LONG_NAME(get_blocksize) (m, n, bl, bu, 0, blk) ;

        dsize = blk[0]*blk[1] > blk[2]*blk[3] ? blk[0]*blk[1] : blk[2]*blk[3] ;
        work = (COV_Entry *) malloc(2 * dsize * sizeof(COV_Entry)) ;
        if (!work)
        {
            free(AllA) ;
            free(mtemp) ;
            free(diag) ;
            printf("Out of memory") ;
            return (0) ;
        }
        iblk = blk ;
    }
    else if (testmem)
    {
        iblk = NULL ;
        work = NULL ;
        my_tries = 0 ;
        /* Make sure the MALLOC fails */
        err = PIRO_BAND(reduce) (iblk, m, n, m, bl, bu, (Entry *) A, ldab, D,
            E+1, (Entry *) U, m, (Entry *) V, n, (Entry *) C, m,
            (Entry *) work, sym) ;
        my_tries = -1 ;
        if (err != PIRO_BAND_OUT_OF_MEMORY)
        {
            printf("Out of Memory test failed "ID" \n", err) ;
            return (0) ;
        }
    }
    else
    {
        /* Override the default blocksize */
        iblk = tblk ;
        work = NULL ;
    }

    /* Reduce the input matrix to the bidiagonal form */
    err = PIRO_BAND(reduce) (iblk, m, n, m, bl, bu, (Entry *) A, ldab, D, E+1,
        (Entry *) U, m, (Entry *) V, n, (Entry *) C, m, (Entry *) work, sym) ;
    if (err != PIRO_BAND_OK)
    {
        free(AllA) ;
        free(mtemp) ;
        free(diag) ;
        if (work != NULL) free(work) ;
        printf("test Band reduction failed "ID"\n", err) ;
        return (0) ;
    }

#if PRINT_TEST_OUTPUT
    PIRO_BAND(print_output) (m, n, D, E, U, V, C) ;
#endif

    /* Check that the A - (U * S * V') is small enough, and check if C=U'  */
    if (!PIRO_BAND(chk_output) (m, n, bu, bl, ldab, sym, D, E, U, V, C, T,0))
    {
        free(AllA) ;
        free(mtemp) ;
        free(diag) ;
        if (work != NULL) free(work) ;
        printf("Chk_output failed\n") ;
        return (0) ;
    }

    /* repeat, but do not provide C */
    memcpy(A, T, n * ldab * sizeof(COV_Entry)) ;
    err = PIRO_BAND(reduce) (iblk, m, n, m, bl, bu, (Entry *) A, ldab, D, E+1,
        (Entry *) U, m, (Entry *) V, n, NULL, m, (Entry *) work, sym);

    /* Check that the A - (U * S * V') is small enough.  Do not check C */
    if (!PIRO_BAND(chk_output) (m, n, bu, bl, ldab, sym, D, E, U, V, NULL,
        T, 0))
    {
        printf ("sym %g\n", (double) sym) ;
        printf("Chk_output (2) failed\n") ;
        abort ( ) ;
        /*
        free(AllA) ;
        free(mtemp) ;
        free(diag) ;
        if (work != NULL) free(work) ;
        return (0) ;
        */
    }

    /* Free Allocated space */
    free(AllA) ;
    free(mtemp) ;
    free(diag) ;
    if (work != NULL) free(work) ;
    return (1) ;
}

/* ========================================================================== */
/* === piro_band_test_error ================================================= */
/* ========================================================================== */

/*
 * Tests all the error cases for the band reduction routine using a matrix
 * stored in a file.
 * Input:
 *  m, n, bl, bu:  matrix is m-by-n with lower and upper bandwidths bl and bu.
 *
 * returns 1 if OK, 0 if failure
 * */
static Int PIRO_BAND(test_error)
(
    Int m,
    Int n,
    Int bl,
    Int bu
)
{
    Int ldab, blk[4], dsize, sym, err, ok = 1 ;
    COV_Entry *A , *work, *U, *V, *C ;
    Entry *D, *E, *diag ;
    Entry *A1 , *work1, *U1, *V1, *C1 ;
    COV_Entry *mtemp ;

    sym = 0 ;

    printf("m = "ID", n= "ID", bl="ID", bu="ID"\n", m, n, bl, bu) ;

    ldab = bl+bu+1 ;

    A = (COV_Entry *) malloc(n * (ldab) * sizeof(COV_Entry)) ;
    if (!A)
    {
        printf("Out of memory") ;
        return (0) ;
    }

    /* Get the numerical values from the file */
    if (!PIRO_BAND(get_random_matrix) (ldab, n, A, sym, bu))
    {
        free (A) ;
        printf("Error in generating matrix\n") ;
        return (0) ;
    }

    /* Check error in get_blocksize */
    err = PIRO_BAND_LONG_NAME(get_blocksize) (-1, n, bl, bu, 0, blk) ;
    ok = ok && (err == PIRO_BAND_M_INVALID) ;

    err = PIRO_BAND_LONG_NAME(get_blocksize) (m, -1, bl, bu, 0, blk) ;
    ok = ok && (err == PIRO_BAND_N_INVALID) ;

    err = PIRO_BAND_LONG_NAME(get_blocksize) (m, n, -1, bu, 0, blk) ;
    ok = ok && (err == PIRO_BAND_BL_INVALID) ;

    err = PIRO_BAND_LONG_NAME(get_blocksize) (m, n, bl, 0, 0, blk) ;
    ok = ok && (err == PIRO_BAND_BU_INVALID) ;

    err = PIRO_BAND_LONG_NAME(get_blocksize) (m, n, bl, bu, 0, NULL) ;
    ok = ok && (err == PIRO_BAND_BLKSIZE_INVALID) ;

    /* Get the suggested block size for the given matrix */
    /* Always use wantuv = 0, as wantuv = 1 will be covered by smaller
     * matrices anyway. */
    PIRO_BAND_LONG_NAME(get_blocksize) (m, n, bl, bu, 0, blk) ;

    /* Allocate space for the singular vectors */
    if (!PIRO_BAND(init_input) (m, n, &mtemp, &U, &V, &C))
    {
        free(A) ;
        printf("Error in initializing input \n") ;
        return (0) ;
    }

    /* Allocate memory for the workspace, bidiagonal (D, E). */
    diag = (Entry *) malloc( 2 * MIN(m, n) * sizeof(Entry)) ;
    if (!diag)
    {
        printf("Out of memory") ;
        free(A) ;
        free(mtemp) ;
        return (0) ;
    }
    D = diag ;
    E = diag + MIN(m, n) ;

    dsize = blk[0]*blk[1] > blk[2]*blk[3] ? blk[0]*blk[1] : blk[2]*blk[3] ;
    work = (COV_Entry *) malloc(2 * dsize * sizeof(COV_Entry)) ;
    if (!work)
    {
        free(A) ;
        free(mtemp) ;
        free(diag) ;
        printf("Out of memory") ;
        return (0) ;
    }

    A1 = (Entry *) A ;
    U1 = (Entry *) U ;
    V1 = (Entry *) V ;
    C1 = (Entry *) C ;
    work1 = (Entry *) work ;

    err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, D, E+1, U1,
        m, V1, n, C1, m, work1, sym) ;
    ok = ok && (err == PIRO_BAND_OK) ;

    /* Test all the error cases */
    E[0] = 0.0 ;
    if (m > 0 && n > 0)
    {
        err = PIRO_BAND(reduce) (blk, -1, n, m, bl, bu, A1, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_M_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, -1, m, bl, bu, A1, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_N_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, -1, bl, bu, A1, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_NRC_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, -1, bu, A1, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_BL_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, -1, A1, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_BU_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, NULL, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_A_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, bl+bu, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_LDAB_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, NULL, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_B1_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, D, NULL, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_B2_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, D, E+1, U1,
            -1, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_LDU_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, D, E+1, U1,
            m, V1, -1, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_LDV_INVALID) ;

        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, D, E+1, U1,
            m, V1, n, C1, -1, work1, sym) ;
        ok = ok && (err == PIRO_BAND_LDC_INVALID) ;

        /* empty matrix: this is not an error condition */
        err = PIRO_BAND(reduce) (blk, 0, 0, m, bl, bu, A1, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_OK) ;

        blk[0] = 0 ;
        blk[1] = 0 ;
        err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, D, E+1, U1,
            m, V1, n, C1, m, work1, sym) ;
        ok = ok && (err == PIRO_BAND_BLKSIZE_INVALID) ;

        if (m != n)
        {
            err = PIRO_BAND(reduce) (blk, m, n, m, bl, bu, A1, ldab, D, E+1, U1,
                m, V1, n, C1, m, work1, 1) ;
            ok = ok && (err == PIRO_BAND_SYM_INVALID) ;
        }
    }

    if (!ok)
    {
        printf ("piro_band_reduce test failed\n") ;
    }

    free(A) ;
    free(mtemp) ;
    free(diag) ;
    free(work) ;
    return (ok) ;
}

/* ========================================================================== */
/* === piro_band_lapack_test ================================================ */
/* ========================================================================== */

/*
 * Test the LAPACK like wrappers for the band reduction routine using a matrix
 * stored in a file.
 * Input:
 *  m, n, bl, bu:  matrix is m-by-n with lower and upper bandwidths bl and bu.
 *  sym      - 0/1. A flag that is one if the input matric is symmetric and
 *             zero if it is unsymmetric.
 *  ul       - A flag to indicate the use of the upper lower triangular part
 *             when the input matrix is symmetric.
 *             1 - upper triangular part. 2 - lower triangular part.
 *  testmem  - 0/1. A flag to turn on the memory allocation tests.
 *
 *  returns 1 if OK, 0 if failure.
 * */
static Int PIRO_BAND(lapack_test)       /* return 1 if OK, 0 if failure */
(
    Int m,
    Int n,
    Int bl,
    Int bu,
    Int sym,            /* create a symmetric test matrix */
    Int ul,
    Int testmem,        /* test memory allocation */
    int *ntests
)
{
    Int i, ldab  ;
    COV_Entry *A , *work, *U, *V, *C, *T, *A2 ;
    Entry *D, *E, *diag ;
    COV_Entry *AllA, *Aptr ;
    COV_Entry *mtemp ;
    Int info, ok ;
    char vect[1] ;
    char uplo[1] ;
    Int sb ;
    Int trial ;

    vect[0] = 'B' ;
    printf("m = "ID", n= "ID", bl="ID", bu="ID"\n", m, n, bl, bu) ;
    (*ntests) ++ ;

    /* Store the upper/lower triangular part of the symmetric matrix */
    if (sym && ul == 1)
    {
        bl = 0 ;
        sb = bu ;
    }
    else if (sym && ul == 2)
    {
        bu = 0 ;
        sb = bl ;
    }

    ldab = bl+bu+1 ;

    AllA = (COV_Entry *) malloc(3 * n * (ldab) * sizeof(COV_Entry)) ;
    if (!AllA)
    {
        printf("Out of memory") ;
        return (0) ;
    }

    /* Allocate memory three times the required memory for the matrix and
     * keep a copy as band reduction will destroy the input matrix. */
    Aptr = AllA ;
    A = Aptr ;
    Aptr += (n * ldab) ;
    T = Aptr ;
    Aptr += (n * ldab) ;
    A2 = Aptr ;

    /* Get the numerical values from the file */
    if (!PIRO_BAND(get_random_matrix) (ldab, n, A, sym, bu))
    {
        free(AllA) ;
        printf("Error in generating matrix \n") ;
        return (0) ;
    }
    memcpy(T, A, n * ldab * sizeof(COV_Entry)) ;
    memcpy(A2, A, n * ldab * sizeof(COV_Entry)) ;


    /* Allocate space for the singular vectors */
    if (!PIRO_BAND(init_input) (m, n, &mtemp, &U, &V, &C))
    {
        free(AllA) ;
        printf("Error in initializing input \n") ;
        return (0) ;
    }

    /* Allocate memory for the workspace, bidiagonal (D, E). */
    diag = (Entry *) malloc( 2 * MIN(m, n) * sizeof(Entry)) ;
    if (!diag)
    {
        printf("Out of memory") ;
        free(A) ;
        free(mtemp) ;
        return (0) ;
    }
    D = diag ;
    E = diag + MIN(m, n) ;

    work = NULL ;
    ok = 1 ;

    if (sym == 0)
    {
        /* Reduce unsymmetric matrix to bidiagonal form */
        E[0] = 0.0 ;
        PIRO_BAND_LAPACK_NAMES(vect, m, n, m, bl, bu, (Entry *) A, ldab,
            D, E+1, (Entry *) U, m, (Entry *) V, n, (Entry *) C,
            m, (Entry *) work, &info) ;
        if (info != 0)
        {
            printf("lapack_test Band reduction failed "ID"\n", info) ;
            ok = 0 ;
        }

#ifdef PRINT_TEST_OUTPUT
        PIRO_BAND(print_output) (m, n, D, E, U, V, C) ;
#endif

        /* Check that the A - (U * S * V') is small enough, and check C=U' */
        if (!PIRO_BAND(chk_output) (m, n, bu, bl, ldab, sym, D, E, U, V, C,
            T, 1))
        {
            free(AllA) ;
            free(mtemp) ;
            free(diag) ;
            printf("Chk_output failed\n") ;
            return (0) ;
        }

        /* Unsymmetric error cases are tested for each matrix. One such
         * test is enough for coverage.
         * */
        E[0] = 0.0 ;
        vect[0] = 'A' ;     /* error case: vect is invalid */
        PIRO_BAND_LAPACK_NAMES(vect, m, n, m, bl, bu, (Entry *) A, ldab,
            D, E+1, (Entry *) U, m, (Entry *) V, n, (Entry *) C,
            m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_VECT_INVALID) ;

        vect[0] = 'B' ;     /* error case: Q is invalid */
        PIRO_BAND_LAPACK_NAMES(vect, m, n, m, bl, bu, (Entry *) A, ldab,
            D, E+1, (Entry *) NULL, m, (Entry *) V, n, (Entry *) C,
            m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_U_INVALID) ;

        /* error case: PT is invalid */
        PIRO_BAND_LAPACK_NAMES(vect, m, n, m, bl, bu, (Entry *) A, ldab,
            D, E+1, (Entry *) U, m, (Entry *) NULL, n, (Entry *) C,
            m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_V_INVALID) ;

        /* error case: C is invalid */
        if (m > 0)
        {
            PIRO_BAND_LAPACK_NAMES(vect, m, n, m, bl, bu, (Entry *) A, ldab,
                D, E+1, (Entry *) U, m, (Entry *) V, n, (Entry *) NULL,
                m, (Entry *) work, &info) ;
            ok = ok && (info == PIRO_BAND_C_INVALID) ;
        }

        /* error case: m is invalid */
        PIRO_BAND_LAPACK_NAMES(vect, -1, n, m, bl, bu, (Entry *) A, ldab,
            D, E+1, (Entry *) U, m, (Entry *) V, n, (Entry *) C,
            m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_M_INVALID) ;

        /* Make sure the MALLOC fails */
        for (trial = 0 ; trial <= 2 ; trial++)
        {
            my_tries = trial ;
            PIRO_BAND_LAPACK_NAMES(vect, m, n, m, bl, bu, (Entry *) A, ldab, D,
                E+1, (Entry *) U, m, (Entry *) V, n,
                (Entry *) C, m, (Entry *) work, &info) ;
            ok = ok &&
                (info == PIRO_BAND_OK || info == PIRO_BAND_OUT_OF_MEMORY) ;
        }
        my_tries = -1 ; /* normal case */
    }
    else
    {
        /* Reduce the symmetric matrix to tridiagonal form */
        vect[0] = 'V' ;
        uplo[0] = (ul == 1 ? 'U' : 'L') ;

        E[0] = 0.0 ;
        PIRO_BAND_LAPACK_SYM_NAMES(vect, uplo, n, sb, (Entry *) A, ldab, D,
            E+1, (Entry *) U, m, (Entry *) work, &info) ;
        if (info != 0)
        {
            printf("lapack_test sym Band reduction failed "ID"\n", info) ;
            ok = 0 ;
        }

#ifdef PRINT_TEST_OUTPUT
        PIRO_BAND(print_output) (m, n, D, E, U, V, C) ;
#endif

        /* Copy U to V and transpose it */
        for ( i = 0 ; i < n*n ; i++)
        {
            V[i] = U[i] ;
        }
        PIRO_BAND(inplace_conjugate_transpose) (n, (Entry *) V, n) ;

        /* Check that the A - (U * S * V') is small enough.  Do not check C */
        if (!PIRO_BAND(chk_output) (m, n, bu, bl, ldab, sym, D, E, U, V, NULL,
            T, 1))
        {
            free(AllA) ;
            free(mtemp) ;
            free(diag) ;
            printf("Chk_output failed\n") ;
            return (0) ;
        }

        /* Test the case where U will not be initialized too */
        vect[0] = 'U' ;
        PIRO_BAND(set_to_eye) (n, n, (Entry *) U) ;
        E[0] = 0.0 ;
        PIRO_BAND_LAPACK_SYM_NAMES(vect, uplo, n, sb, (Entry *) A2, ldab, D,
            E+1, (Entry *) U, m, (Entry *) work, &info) ;
        if (info != 0)
        {
            printf("lapack_test sym Band reduction failed "ID"\n", info) ;
            ok = 0 ;
        }

#ifdef PRINT_TEST_OUTPUT
        PIRO_BAND(print_output) (m, n, D, E, U, V, C) ;
#endif

        /* Copy U to V and transpose it */
        for ( i = 0 ; i < n*n ; i++)
        {
            V[i] = U[i] ;
        }
        PIRO_BAND(inplace_conjugate_transpose) (n, (Entry *) V, n) ;

        /* Check that the A - (U * S * V') is small enough.  Do not check C */
        if (!PIRO_BAND(chk_output) (m, n, bu, bl, ldab, sym, D, E, U, V, NULL,
            T, 1))
        {
            free(AllA) ;
            free(mtemp) ;
            free(diag) ;
            printf("Chk_output failed\n") ;
            return (0) ;
        }

        /* Symmetric error cases are tested for each matrix. One such
         * test is enogh for coverage.
         * */
        E[0] = 0.0 ;

        vect[0] = 'V' ; /* error case : Q is NULL */
        PIRO_BAND_LAPACK_SYM_NAMES(vect, uplo, n, sb, (Entry *) A, ldab, D,
            E+1, (Entry *) NULL, m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_U_INVALID) ;

        vect[0] = 'A' ; /* error case: vect is invalid */
        PIRO_BAND_LAPACK_SYM_NAMES(vect, uplo, n, sb, (Entry *) A, ldab, D,
            E+1, (Entry *) U, m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_VECT_INVALID) ;

        vect[0] = 'V' ;
        uplo[0] = 'A' ; /* error case: uplo is invalid */
        PIRO_BAND_LAPACK_SYM_NAMES(vect, uplo, n, sb, (Entry *) A, ldab, D,
            E+1, (Entry *) U, m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_UPLO_INVALID) ;

        /* error case */
        uplo[0] = (ul == 1 ? 'U' : 'L') ;
        PIRO_BAND_LAPACK_SYM_NAMES(vect, uplo, n, -1, (Entry *) A, ldab, D,
            E+1, (Entry *) U, m, (Entry *) work, &info) ;
        ok = ok && (info == PIRO_BAND_BU_INVALID) ;

        /* Make sure the MALLOC fails */
        vect[0] = 'V' ;
        for (trial = 0 ; testmem && trial <= 2 ; trial++)
        {
            my_tries = trial ;
            PIRO_BAND_LAPACK_SYM_NAMES(vect, uplo, n, sb, (Entry *) A, ldab, D,
                                    E+1, (Entry *) U, m, (Entry *) work, &info);
            ok = ok &&
                (info == PIRO_BAND_OK || info == PIRO_BAND_OUT_OF_MEMORY) ;
        }
        my_tries = -1 ; /* normal case */
    }

    free(AllA) ;
    free(mtemp) ;
    free(diag) ;
    if (!ok) printf ("LAPACK test failure, info %g\n", (double) info) ;
    return (ok) ;
}

/* ========================================================================== */
/* === piro_band_call_givens ================================================ */
/* ========================================================================== */

/* Call the givens rotations function to cover special cases and compare the
 * result d returned to the fist component after applying the rotations. For
 * the PIROBAND_FLOAT version we can compare the result to the result obtained
 * by calling the double version of the code. But that part is commented out
 * now.
 */
static void PIRO_BAND(call_givens)
(
    Entry *f,
    Entry *g,
    PIRO_BAND_Common *Common
)
{
    Entry c[2], s[2], r[2] ;
    Entry temp[2] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[2] ;
#else
    Entry temp2[2] ;
#endif

    /* Find the givens rotations */
    PIRO_BAND(givens) (f, g, c, s, r, Common) ;
    APPLY_GIVENS(f, g, c, s, temp, conjs) ;

    /* Compare the r from the Givens and the result of applying the givens
     * rotations to f and g */
#ifdef PIROBAND_COMPLEX
    printf("r - (G * [f; g]) (1)=%g %g\n", r[0]-f[0], r[1]-f[1] ) ;
#else
    printf("r - (G * [f; g]) (1)= %g\n", r[0]-f[0] ) ;
#endif

}

/* ========================================================================== */
/* === piro_band_test_givens ================================================ */
/* ========================================================================== */

/* Test the Givens rotation functions for special cases */
static void PIRO_BAND(test_givens) ( void )
{
    Entry safemin, safemin2, safemx2, safemx ;
    Entry eps ;
    Int esfmn2 ;
    Entry f[2], g[2] ;
    PIRO_BAND_Common Common ;

    /* Set simpler min and max for coverage */
    safemin2 = sqrt(PIRO_BAND_Entry_MIN/PIRO_BAND_EPS) ;
    safemx2 = 1/safemin2 ;
    safemin = PIRO_BAND_Entry_MIN ;
    safemx = COV_DBL_MAX ;

    /* Set the safe minimum and safe maximum for the givens */
    eps = PIRO_BAND_EPS/2 ;
    /* logf and powf are C99 requirements. Use log and pow always */
    esfmn2  = log(safemin/eps) / log(2) / 2 ;
    Common.safemin2 = pow(2, esfmn2) ;
    Common.safemx2 = 1/Common.safemin2 ;

    f[0] = safemx2 * 1.1 ;
    f[1] = 3.0 ;
    g[0] = 1009.0 ;
    g[1] = 4.0 ;
    PIRO_BAND(call_givens) (f, g, &Common) ;

    f[0] = 0.0 ;
    f[1] = 0.0 ;
    g[0] = 1009.0 ;
    g[1] = 133.0 ;
    PIRO_BAND(call_givens) (f, g, &Common) ;

    f[0] = safemin2 * 0.9 ;
    f[1] = safemin2 * 0.1 ;
    g[0] = safemin2 * 0.8 ;
    g[1] = safemin2 * 0.1 ;
    PIRO_BAND(call_givens) (f, g, &Common) ;

    f[0] = safemin2 * 1.1 ;
    f[1] = 0.0 ;
    g[0] = 9 * (safemin2/safemin) ;
    g[1] = 0 ;
    PIRO_BAND(call_givens) (f, g, &Common) ;

    f[0] = 1.1 ;
    f[1] = 0.0 ;
    g[0] = sqrt(safemx) ;
    g[1] = 0.0 ;
    PIRO_BAND(call_givens) (f, g, &Common) ;
}

/* ========================================================================== */
/* === piro_band_run_tests ================================================== */
/* ========================================================================== */

/* Test band reduction for various input matrices for all the matrices */
Int PIRO_BAND(run_tests) ( void )    /* returns 1 if OK, 0 if failure */
{
    Int blks[4], n, b ;
    int pass = 0, ntests = 0 ;

    pass += PIRO_BAND(test) (10, 10, 3, 4, 0, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (10, 10, 4, 4, 0, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (300, 300, 101, 110, 0, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (100, 100, 40, 40, 1, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (10, 20, 7, 8, 0, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (10, 20, 8, 8, 0, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (20, 20, 1, 2, 0, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (10, 21, 5, 5, 0, 0, NULL, &ntests) ;
    pass += PIRO_BAND(test) (15, 15, 0, 1, 0, 0, NULL, &ntests) ;

    blks[0] = 2 ;
    blks[1] = 3 ;
    blks[2] = 2 ;
    blks[3] = 3 ;
    pass += PIRO_BAND(test) (151, 150, 90, 100, 0, 0, blks, &ntests) ;

    blks[0] = 2 ;
    blks[1] = 3 ;
    blks[2] = 3 ;
    blks[3] = 1 ;
    pass += PIRO_BAND(test) (10, 21, 5, 5, 0, 0, blks, &ntests) ;

    pass += PIRO_BAND(test) (10, 21, 5, 5, 0, 1, NULL, &ntests) ;

    PIRO_BAND(test_error) (100, 100, 40, 40) ;
    PIRO_BAND(test_error) (10, 20, 7, 8) ;

    printf("--------- Testing LAPACK style interface -------------\n") ;
    pass += PIRO_BAND(lapack_test) (10, 10, 3, 4, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (10, 10, 4, 4, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (300, 300, 101, 110, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (150, 150, 100, 90, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (150, 150, 90, 100, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (20, 20, 7, 7, 1, 1, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (20, 20, 7, 7, 1, 2, 1, &ntests) ;
    pass += PIRO_BAND(lapack_test) (10, 10, 5, 0, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (10, 10, 1, 0, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (10, 10, 0, 0, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (10, 10, 0, 5, 0, 0, 0, &ntests) ;
    pass += PIRO_BAND(lapack_test) (10, 10, 0, 0, 1, 1, 1, &ntests) ;

    printf("---------- Testing Givens rotations -------------- \n") ;
    PIRO_BAND(test_givens) () ;

    /* test blocksizes for large matrices */
    n = 10000 ;
    b = 1000 ;
    pass += (PIRO_BAND_LONG_NAME(get_blocksize) (n, n, b, b, 0, blks) ==
        PIRO_BAND_OK) ;
    ntests ++ ;

    if (pass != ntests)
    {
        printf ("Test failure: only %g of %g tests passed\n",
            (double) pass, (double) ntests) ;
        return (0) ;
    }

    return (1) ;
}
