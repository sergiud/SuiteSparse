/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_main.c ==================================== */
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
 * One main file that includes all 8 possible versions of piro_band methods.
 * The methods that are dependent only on PIROBAND_LONG are not included here.
 * They are compiled twice.
 */

/* ------------------------ Double  ---------------------------------------- */
#include "piro_band_internal.h"
#include "piro_band.c" /* dri */
#include "piro_band_givens.c"

#define PIROBAND_COMPLEX      /* dci */
#include "piro_band_internal.h"
#include "piro_band.c"
#include "piro_band_givens.c"
#undef PIROBAND_COMPLEX
/*It is enough to include piro_band_blocksize.h and piro_band_uv_update.h
 * twice once each for integer and long versions of the code
 * */
#undef PIRO_BAND_UV_UPDATE_H
/* Introduce a new typedef for uv_update by undefining L_UV for the long case
 * of U and V update.
 * */
#undef L_UV

#define PIROBAND_LONG         /* drl */
#include "piro_band_internal.h"
#include "piro_band.c"
#include "piro_band_givens.c"

#define PIROBAND_COMPLEX      /* dcl */
#include "piro_band_internal.h"
#include "piro_band.c"
#include "piro_band_givens.c"
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG

/* ------------------------ Single  ---------------------------------------- */
#define PIROBAND_FLOAT         /* sri */
#include "piro_band_internal.h"
#include "piro_band.c"
#include "piro_band_givens.c"

#define PIROBAND_COMPLEX       /* sci */
#include "piro_band_internal.h"
#include "piro_band.c"
#include "piro_band_givens.c"
#undef PIROBAND_COMPLEX

#define PIROBAND_LONG          /* srl */
#include "piro_band_internal.h"
#include "piro_band.c"
#include "piro_band_givens.c"

#define PIROBAND_COMPLEX       /* scl */
#include "piro_band_internal.h"
#include "piro_band.c"
#include "piro_band_givens.c"
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG
#undef PIROBAND_FLOAT
