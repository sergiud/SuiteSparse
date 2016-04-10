/* ========================================================================== */
/* === PIRO_BAND/Tcov/memory.c ============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Memory allocation testing for piro_band functions. my_malloc2, pretend to
 * fail if my_tries goes to zero.  No failure occurs if my_tries is negative.
 */

#include <stdlib.h>
#include "piro_band_internal.h"
#include "mem.h"

/* ========================================================================== */
/* === my_tries ============================================================= */
/* ========================================================================== */

extern int my_tries ;

int my_tries = -1 ; /* a global variable */

/* ========================================================================== */
/* === my_malloc2 =========================================================== */
/* ========================================================================== */

void *my_malloc2 (size_t size)
{
    void *p ;
    if (my_tries == 0)
    {
        /* pretend to fail */
        return (NULL) ;
    }
    if (my_tries > 0)
    {
        my_tries-- ;
    }
    p = malloc (size) ;
    return (p) ;
}

/* ========================================================================== */
/* === my_free2 ============================================================= */
/* ========================================================================== */

void my_free2 (void *p)
{
    free (p) ;
}

