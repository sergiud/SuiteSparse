/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_blocksize_main.c ========================== */
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
 * Includes 2 possible versions of piro_band_blocksize functions.
 */

/* ------------------------ Double  ---------------------------------------- */
#include "piro_band_internal.h"
#include "piro_band_blocksize.c"
#undef PIRO_BAND_BLOCKSIZE_H
#undef L_UV

#define PIROBAND_LONG
#include "piro_band_internal.h"
#include "piro_band_blocksize.c"
#undef PIROBAND_LONG
#undef PIRO_BAND_BLOCKSIZE_H

