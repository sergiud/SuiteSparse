/* ========================================================================== */
/* === PIRO_BAND/Include/piro_band_blocksize.h ============================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Function to determine the recommended blocksize for a given banded matrix of
 * size mxn, a lower bandwidth of bl and upper bandwidth of bu.  Output is in
 * blk where blk[0] and blk[1]  is the #columns and #rows in the block for the
 * reducing the upper bandwidth. blk[2] and blk[2] is the #rows and #columns in
 * the block for reducing the lower bandwidth.
 *
 * This file should not be #include'd in user programs.
 * */
#ifndef PIRO_BAND_BLOCKSIZE_H
#define PIRO_BAND_BLOCKSIZE_H

int PIRO_BAND_LONG_NAME(get_blocksize)
(
    Int m,              /* #rows in the input matrix */
    Int n,              /* #columns in the input matrix */
    Int bl,             /* lower bandwidth of the input matrix */
    Int bu,             /* upper bandwidth of the input matrix */
    Int wantuv,         /* flag : 0/1 whether U or V is needed */
    Int *blk            /* Output, recommended block size. Shoule be atleast of
                         * size 4. */
) ;

#endif /* PIRO_BAND_BLOCKSIZE_H */
