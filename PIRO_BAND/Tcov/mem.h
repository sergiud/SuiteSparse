/* ========================================================================== */
/* === PIRO_BAND/Tcov/mem.h ================================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Functions for overriding malloc and free and testing memory allocation in
 * code coverage
 * */

void *my_malloc2 (size_t size) ;
void my_free2 (void *p) ;
