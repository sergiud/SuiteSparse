/* ========================================================================== */
/* === PIRO_BAND/MATLAB/piro_band_qr_main.c ================================= */
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
 * Includes 2 possible versions of piro_band_qr methods.  The other 6 methods
 * are not required for MATLAB interfaces.
 */

/* ------------------------ Double  ---------------------------------------- */
#define PIROBAND_LONG         /* drl */
#include "piro_band_matlab.h"
#include "piro_band_qr.c"
#undef PIRO_BAND_MATLAB_H

#define PIROBAND_COMPLEX      /* dcl */
#include "piro_band_matlab.h"
#include "piro_band_qr.c"
#undef PIROBAND_COMPLEX
#undef PIROBAND_LONG

