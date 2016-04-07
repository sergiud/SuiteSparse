/* ========================================================================== */
/* === PIRO_BAND/Tcov/piro_band_cover.h ===================================== */
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
 * Assign the correct names for functions that allow LAPACK style interface
 * for bidiagonal and tridiagonal band reduction of band matrice based on
 * whether PIROBAND_LONG, PIROBAND_FLOAT and PIROBAND_COMPLEX are defined or
 * not.
 */

#ifndef PIRO_BAND_COVER_H
#define PIRO_BAND_COVER_H

#undef PIRO_BAND_LAPACK_NAMES
#undef PIRO_BAND_LAPACK_SYM_NAMES

#ifdef PIROBAND_LONG

#ifdef PIROBAND_FLOAT

#ifdef PIROBAND_COMPLEX

#define PIRO_BAND_LAPACK_NAMES piro_band_cgbbrd_l
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_chbtrd_l

#else /* PIROBAND_COMPLEX */

#define PIRO_BAND_LAPACK_NAMES piro_band_sgbbrd_l
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_ssbtrd_l

#endif /* PIROBAND_COMPLEX */

#else /* PIROBAND_FLOAT */

#ifdef PIROBAND_COMPLEX

#define PIRO_BAND_LAPACK_NAMES piro_band_zgbbrd_l
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_zhbtrd_l

#else /* PIROBAND_COMPLEX */

#define PIRO_BAND_LAPACK_NAMES piro_band_dgbbrd_l
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_dsbtrd_l

#endif /* PIROBAND_COMPLEX */

#endif /* PIROBAND_FLOAT */

#else /* PIROBAND_LONG */

#ifdef PIROBAND_FLOAT

#ifdef PIROBAND_COMPLEX

#define PIRO_BAND_LAPACK_NAMES piro_band_cgbbrd
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_chbtrd

#else /* PIROBAND_COMPLEX */

#define PIRO_BAND_LAPACK_NAMES piro_band_sgbbrd
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_ssbtrd

#endif /* PIROBAND_COMPLEX */

#else /* PIROBAND_FLOAT */

#ifdef PIROBAND_COMPLEX

#define PIRO_BAND_LAPACK_NAMES piro_band_zgbbrd
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_zhbtrd

#else /* PIROBAND_COMPLEX */

#define PIRO_BAND_LAPACK_NAMES piro_band_dgbbrd
#define PIRO_BAND_LAPACK_SYM_NAMES piro_band_dsbtrd

#endif /* PIROBAND_COMPLEX */

#endif /* PIROBAND_FLOAT */

#endif /* PIROBAND_LONG */

#endif /* PIRO_BAND_COVER_H */
