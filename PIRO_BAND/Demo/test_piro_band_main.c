/* ========================================================================== */
/* === PIRO_BAND/Test/test_piro_band_main.c ================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Main Source file for testing piro_band_reduce . This just includes the other
 * source files multiple times for the various versions we support.
 *  */

#include "piro_band.h"
#include "piro_band_lapack.h"
#include "piro_band_lapack_internal.h"

/* ------------------------ Double  ---------------------------------------- */

/* dri: double, real, int */
#include "piro_band_internal.h"
#include "test_piro_band.c" 
#include "piro_band_testutils.c"

/* dci: double, complex, int */
#define PIROBAND_COMPLEX
#include "piro_band_internal.h"
#include "test_piro_band.c"
#include "piro_band_testutils.c"
#undef PIROBAND_COMPLEX
#undef PIRO_BAND_BLOCKSIZE_H
#undef PIRO_BAND_UV_UPDATE_H
#undef L_UV

/* drl: double, real, long */
#define PIROBAND_LONG
#include "piro_band_internal.h"
#include "test_piro_band.c"
#include "piro_band_testutils.c"

/* dcl: double, complex, long */
#define PIROBAND_COMPLEX
#include "piro_band_internal.h"
#include "test_piro_band.c"
#include "piro_band_testutils.c"
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG

/* ------------------------ Single  ---------------------------------------- */

/* sri: float, real, int */
#define PIROBAND_FLOAT
#include "piro_band_internal.h"
#include "test_piro_band.c"
#include "piro_band_testutils.c"

/* sci: float, complex, int */
#define PIROBAND_COMPLEX
#include "piro_band_internal.h"
#include "test_piro_band.c"
#include "piro_band_testutils.c"
#undef PIROBAND_COMPLEX

/* srl: float, real, long */
#define PIROBAND_LONG
#include "piro_band_internal.h"
#include "test_piro_band.c"
#include "piro_band_testutils.c"

/* scl: float, complex, long */
#define PIROBAND_COMPLEX
#include "piro_band_internal.h"
#include "test_piro_band.c"
#include "piro_band_testutils.c"
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG
#undef PIROBAND_FLOAT

/* Call all the tests explicitly */
int main(void)
{
    int pass = 0 ;
    printf("------------- Tests for piro_band_reduce ------- \n") ;
    printf("------------- Tests for double/real/integer ------- \n") ;
    pass += piro_band_run_tests_dri() ;
    printf("------------- Tests for double/complex/integer ------- \n") ;
    pass += piro_band_run_tests_dci() ;
    printf("------------- Tests for double/real/long ------- \n") ;
    pass += piro_band_run_tests_drl() ;
    printf("------------- Tests for double/complex/long ------- \n") ;
    pass += piro_band_run_tests_dcl() ;
    printf("------------- Tests for float/real/integer ------- \n") ;
    pass += piro_band_run_tests_sri() ;
    printf("------------- Tests for float/complex/integer ------- \n") ;
    pass += piro_band_run_tests_sci() ;
    printf("------------- Tests for float/real/long ------- \n") ;
    pass += piro_band_run_tests_srl() ;
    printf("------------- Tests for float/complex/long ------- \n") ;
    pass += piro_band_run_tests_scl() ;
    if (pass == 8)
    {
        printf("PIRO_BAND:  all tests passed\n") ;
    }
    else
    {
        printf("PIRO_BAND:  TEST FAILURE: only %d of 8 tests passed\n", pass) ;
    }
    return 0 ;
}
