/* ========================================================================== */
/* === PIRO_BAND/Test/test_rand.c =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

#include "test_rand.h"

static unsigned long piro_band_next = 1;

/* RAND_MAX assumed to be 32767 */
int my_rand(void)
{
    piro_band_next = piro_band_next * 1103515245 + 12345 ;
    return ((unsigned) (piro_band_next / 65536) % 32768) ;
}

void mysrand(unsigned seed)
{
    piro_band_next = seed ;
}

double xrand (double range)
{
    return ((range * (double) (my_rand ( ))) / MY_RAND_MAX) ;
}


