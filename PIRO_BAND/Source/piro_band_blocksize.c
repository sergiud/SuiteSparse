/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_blocksize.c =============================== */
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
 * size mxn, a lower bandwidth of bl and upper bandwidth of bu.  blk is an
 * array of 4 integers.  Output is in blk where * blk[0] and blk[1]  is the
 * #columns and #rows in the block for the reducing the upper bandwidth. blk[2]
 * and blk[2] is the #rows and #columns in the block for reducing the lower
 * bandwidth.
 *
 * Returns 0 on success, < 0 on failure.  See piro_band.h for error codes.
 * */

#include "piro_band.h"

int PIRO_BAND_LONG_NAME(get_blocksize)
(
    Int m,              /* #rows in the input matrix */
    Int n,              /* #columns in the input matrix */
    Int bl,             /* lower bandwidth of the input matrix */
    Int bu,             /* upper bandwidth of the input matrix */
    Int wantuv,         /* flag : 0/1 whether U or V is needed */
    Int *blk            /* Output, recommended block size. Array of size 4 */
)
{
    double minmn, d_bw, fl ;

    if (m < 0)
    {
        return (PIRO_BAND_M_INVALID) ;
    }
    if (n < 0)
    {
        return (PIRO_BAND_N_INVALID) ;
    }
    if (bl < 0)
    {
        return (PIRO_BAND_BL_INVALID) ;
    }
    if (bu < 1)
    {
        return (PIRO_BAND_BU_INVALID) ;
    }
    if (blk == NULL)
    {
        return (PIRO_BAND_BLKSIZE_INVALID) ;
    }

    blk[0] = 0 ;
    blk[1] = 0 ;
    blk[2] = 0 ;
    blk[3] = 0 ;

    d_bw = bl + bu ;
    minmn = (m <= n ? m : n) ;
    /* Rough estimate of the flops cutoff point that is to choose larger
     * block size. */
    fl = 6 * d_bw * minmn * minmn ;

    /* Find the block size for the upper band */
    if (bu > 1)
    {
        if (wantuv || bu < 16)
        {
            /* Entire row is the block */
            blk[0] = bu-1 ;
            blk[1] = 1 ;
        }
        else
        {
            if (fl > 1e+10 && bu >= 64)
            {
                /* flops high. Use 32x32. We can use 64x64 in very high cases.
                 * But 32x32 is not very bad in those cases too.
                 * */
                blk[0] = 32 ;
                blk[1] = 32 ;
            }
            else
            {
                /* If flops is less than cutoff or bu is small to use 32, use
                 * 8x8. We can use 16x16 in some cases. But 8x8 is not very bad.
                 * We can use higher cutoff for unsymmetric but the current
                 * cutoff is not bad in those cases too.
                 * */
                blk[0] = 8 ;
                blk[1] = 8 ;
            }
        }
    }
    else
    {
        blk[0] = 0 ;
        blk[1] = 0 ;
    }

    /* Find the block size for the lower band */
    if (bl > 1)
    {
        if (wantuv || bl < 16)
        {
            /* Entire row is the block */
            blk[2] = bl ;
            blk[3] = 1 ;
        }
        else
        {
            if (fl > 1e+10 && bl >= 64)
            {
                /* flops high. Use 32x32. We can use 64x64 in very high cases.
                 * But 32x32 is not very bad in those cases too.
                 * */
                blk[2] = 32 ;
                blk[3] = 32 ;
            }
            else
            {
                /* If flops is less than cutoff or bl is small to use 32, use
                 * 8x8. We can use 16x16 in some cases. But 8x8 is not very bad.
                 * We can use higher cutoff for unsymmetric but the current
                 * cutoff is not bad in those cases too.
                 * */
                blk[2] = 8 ;
                blk[3] = 8 ;
            }
        }
    }
    else
    {
        blk[2] = bl ;
        blk[3] = bl ;
    }

    return (PIRO_BAND_OK) ;
}
