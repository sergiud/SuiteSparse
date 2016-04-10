/* ========================================================================== */
/* === PIRO_BAND/Include/piro_band_memory.h ================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* This file should not be #include'd in user programs. */

/* Functions for overriding malloc and free and testing memory allocation in
 * code coverage.

 TODO delete this file
 * */
#ifndef PIRO_BAND_MEMORY_H
#define PIRO_BAND_MEMORY_H
#include <stdlib.h>

#ifdef TCOV

void *my_malloc2 (size_t size) ;
void my_free2 (void *p) ;

#endif
#endif
