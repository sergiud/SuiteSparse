/* ========================================================================== */
/* === PIRO_BAND/Tcov/piro_band_cover_main.c ================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Main source file for test coverage. This just includes the other source files
 * multiple times for the various versions we support.
 *  */

/* Global variable to control the tests for memory allocation */
extern int my_tries ;


/* ------------------------ Double  ---------------------------------------- */

#include "piro_band_internal.h"
#include "piro_band_cover.c"      /* dri */
#include "piro_band_testutils.c" /* dri */
#undef PIRO_BAND_COVER_H

#define PIROBAND_COMPLEX      /* dci */
#include "piro_band_internal.h"
#include "piro_band_cover.c"
#include "piro_band_testutils.c"
#undef PIROBAND_COMPLEX
#undef PIRO_BAND_COVER_H
#undef PIRO_BAND_BLOCKSIZE_H
#undef PIRO_BAND_UV_UPDATE_H
#undef L_UV

#define PIROBAND_LONG         /* drl */
#include "piro_band_internal.h"
#include "piro_band_cover.c"
#include "piro_band_testutils.c"
#undef PIRO_BAND_COVER_H

#define PIROBAND_COMPLEX      /* dcl */
#include "piro_band_internal.h"
#include "piro_band_cover.c"
#include "piro_band_testutils.c"
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG
#undef PIRO_BAND_COVER_H

/* ------------------------ Single  ---------------------------------------- */
#define PIROBAND_FLOAT         /* sri */
#include "piro_band_internal.h"
#include "piro_band_cover.c"
#include "piro_band_testutils.c"
#undef PIRO_BAND_COVER_H

#define PIROBAND_COMPLEX       /* sci */
#include "piro_band_internal.h"
#include "piro_band_cover.c"
#include "piro_band_testutils.c"
#undef PIROBAND_COMPLEX
#undef PIRO_BAND_COVER_H

#define PIROBAND_LONG          /* srl */
#include "piro_band_internal.h"
#include "piro_band_cover.c"
#include "piro_band_testutils.c"
#undef PIRO_BAND_COVER_H

#define PIROBAND_COMPLEX       /* scl */
#include "piro_band_internal.h"
#include "piro_band_cover.c"
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
        printf ("PIRO_BAND:  all tests passed\n") ;
    }
    else
    {
        printf("PIRO_BAND:  TEST FAILURE: only %d of 8 tests passed\n", pass) ;
    }
    return 0 ;
}
