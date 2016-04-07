/* ========================================================================== */
/* === PIRO_BAND/Include/piro_band_uv_update.h ============================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Methods to set the data structure for the update of the right and left
 * singular/eigen vectors.
 * */

/* This file should not be #include'd in user programs. */

#ifndef PIRO_BAND_UV_UPDATE_H
#define PIRO_BAND_UV_UPDATE_H


/* Data Structure for Accumulating Plane rotations in the right and left side
 * */
struct PIRO_BAND_LONG_NAME(uv_update_struct)
{
    Int u_srow_min ;            /* starting row of U for k = k - u_srow_min   */
    Int vt_srow_min ;           /* starting row of VT for k = k - vt_srow_min */

    Int blksize ;               /* block size for fill                        */
    Int t_blksize ;             /* one of the block size for fill when nr > 1 */

    Int vt_lrow_upper ;         /* last row of VT in upper band */
    Int vt_lrow_lower ;         /* last row of VT in lower band */
    Int u_lrow_upper ;          /* last row of U in upper band  */
    Int u_lrow_lower ;          /* last row of U in lower band  */

    Int bw ;                    /* Original bandwidth of A      */
    Int c_lrow ;                /* current last row             */
    Int k ;                     /* current row/column           */
    Int first ;                 /* flag for first iteration     */
} ;

/* Initialize the value for updating U and VT */
void PIRO_BAND_LONG_NAME(init_uv_update)
(
    Int m,                                  /* # rows in A                   */
    Int n,                                  /* #columns in A                 */
    Int bl,                                 /* lower bandwidth of A          */
    Int bu,                                 /* upper bandwidth of A          */
    Int nr,                                 /* # rows in the block           */
    PIRO_BAND_LONG_NAME(uv_update) *upd       /* data structure for the update */
) ;


/* Set the values for updating U and VT after an iteration */
void PIRO_BAND_LONG_NAME(set_uv_update)
(
    Int m,                                  /* # rows in A                   */
    Int n,                                  /* #columns in A                 */
    Int bl,                                 /* lower bandwidth of A          */
    Int bu,                                 /* upper bandwidth of A          */
    Int nr,                                 /* # rows in the block           */
    PIRO_BAND_LONG_NAME(uv_update) *upd     /* data structure for the update */
) ;
#endif
