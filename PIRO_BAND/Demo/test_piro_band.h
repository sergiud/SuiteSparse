/* ========================================================================== */
/* === PIRO_BAND/Test/test_piro_band.h ====================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* All the definitions are based on whether PIROBAND_FLOAT, PIROBAND_COMPLEX
 * and PIRBAND_LONG are defined
 * or not. (except MIN, MAX, and INDEX)
 * C99 style complex is assumed here.
 * */

#undef COV_Entry
#undef COV_CONJ
#undef COV_ABS
#undef Complex
#undef INDEX
#undef MIN
#undef MAX
#undef PRINT_QRVALUES
#undef COV_DBL_MAX
#undef PRINT_TEST_OUTPUT

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Macros that are not dependent on PIRO_BAND_LONG, PIROBAND_FLOAT,
 * PIROBAND_COMPLEX. They will be the same for all 8 versions.
 */
#define INDEX(row, col) ((row)-(col)+bu+((ldab)*(col)))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/* Uncomment/Comment the following line to turn on/off debugging for the test
 * coverage
 * */
/*
#define PRINT_TEST_OUTPUT
*/


/* ========================================================================= */
/* ======= Definitions for complex/real and double/float usage  ============ */
/* ========================================================================= */
#ifdef PIROBAND_FLOAT

#define COV_DBL_MAX FLT_MAX

#ifdef PIROBAND_COMPLEX

#define COV_ABS cabsf
#define COV_CONJ conjf

#else /* PIROBAND_COMPLEX */

#define COV_ABS fabs
#define COV_CONJ

#endif /* PIROBAND_COMPLEX */

#else /* PIROBAND_FLOAT */

#define COV_DBL_MAX DBL_MAX

#ifdef PIROBAND_COMPLEX

#define COV_ABS cabs
#define COV_CONJ conj

#else /* PIROBAND_COMPLEX */

#define COV_ABS fabs
#define COV_CONJ

#endif /* PIROBAND_COMPLEX */

#endif /* PIROBAND_FLOAT */

/* ========================================================================= */
/* ======= Definitions for complex/real usage  ============================= */
/* ========================================================================= */
#ifdef PIROBAND_COMPLEX

#define Complex _Complex

#define PRINT_QRVALUES(str, i, val, i2) printf(str "%0.4f %0.4f \n", i, \
                                        val[i2+0], val[i2+1])

#else /* PIROBAND_COMPLEX */

#define Complex

#define PRINT_QRVALUES(str, i, val, i2) printf(str "%0.4f \n", i, val[i2+0])


#endif /* PIROBAND_COMPLEX */

#define COV_Entry Complex Double


/* ========================================================================= */
/* ======= Utility functions for tests  ==================================== */
/* ========================================================================= */
void PIRO_BAND(svdmult)
(
    COV_Entry *U,
    COV_Entry *VT,
    Entry *B1,
    Entry *B2,
    Int m,
    Int n,
    COV_Entry *A1,
    Int sym
) ;

Entry PIRO_BAND(find_norm)
(
    COV_Entry *A,
    Int m,
    Int n
) ;

Entry PIRO_BAND(find_band_norm)
(
    COV_Entry *A,
    Int ldab,
    Int m,
    Int n,
    Int bl,
    Int bu
) ;

Int PIRO_BAND(get_matrix)
(
    char *pfile,
    Int ldab,
    Int n,
    COV_Entry *A
) ;

Int PIRO_BAND(get_random_matrix)
(
    Int ldab,
    Int n,
    COV_Entry *A,
    Int sym,
    Int bu
) ;

Int PIRO_BAND(init_input)
(
    Int m,
    Int n,
    COV_Entry **temp1,
    COV_Entry **U,
    COV_Entry **V,
    COV_Entry **C
) ;

void PIRO_BAND(identity_matrix)
(
    Int m,
    COV_Entry *C
) ;

Int PIRO_BAND(chk_output)
(
    Int m,
    Int n,
    Int bu,
    Int bl,
    Int ldab,
    Int sym,
    Entry *D,
    Entry *E,
    COV_Entry *U,
    COV_Entry *V,
    COV_Entry *C,
    COV_Entry *A,
    Int lapack
) ;

Int PIRO_BAND(run_tests)
(void ) ;

void PIRO_BAND(print_output)
(
    Int m,
    Int n,
    Entry *D,
    Entry *E,
    COV_Entry *U,
    COV_Entry *V,
    COV_Entry *C
) ;

