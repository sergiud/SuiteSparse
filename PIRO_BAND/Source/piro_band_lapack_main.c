/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_lapack_main.c ============================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/*
 * Includes all 8 possible versions of piro_band_lapack methods. All the
 * functions depend only on the macro PIROBAND_LONG and not on
 * PIROBAND_COMPLEX or PIROBAND_FLOAT.
 * piro_band_lapack.c is included once for the integer and once for long case.
 * But all the 8 possible prototypes for the functions are included so they can
 * be used.
 */

/* ------------------------ Double  ---------------------------------------- */
#include "piro_band_internal.h"

#define PIROBAND_COMPLEX      /* dci */
#include "piro_band_internal.h"
#undef PIROBAND_COMPLEX
/* Introduce a new typedef for uv_update by undefining L_UV for the long case
 * of U and V update
 * */
#undef L_UV

#define PIROBAND_LONG         /* drl */
#include "piro_band_internal.h"

#define PIROBAND_COMPLEX      /* dcl */
#include "piro_band_internal.h"
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG

/* ------------------------ Single  ---------------------------------------- */
#define PIROBAND_FLOAT         /* sri */
#include "piro_band_internal.h"

#define PIROBAND_COMPLEX       /* sci */
#include "piro_band_internal.h"
#include "piro_band_lapack.c" /* int lapack */
#undef PIROBAND_COMPLEX

#define PIROBAND_LONG          /* srl */
#include "piro_band_internal.h"
#undef PIRO_BAND_BLOCKSIZE_H

#define PIROBAND_COMPLEX       /* scl */
#include "piro_band_internal.h"
#include "piro_band_lapack.c" /* long lapack */
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG
#undef PIROBAND_FLOAT
