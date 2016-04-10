/* ========================================================================== */
/* === PIRO_BAND/Test/test_rand.h =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

#ifndef RAND_H
#define RAND_H

#define MY_RAND_MAX 32767

/* RAND_MAX assumed to be 32767 */
int my_rand (void) ;

void mysrand (unsigned seed) ;

double xrand (double range) ;

#endif
