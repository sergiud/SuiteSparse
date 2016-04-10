/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_uv_update.c =============================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* #include "piro_band_internal.h" */
#include "piro_band_uv_update.h"

/* Functions to initialize and set the uv_update data structure for update of U
 * and V when using the blocked algorithm for the equal bandwidth case.
 */

/* ========================================================================== */
/* init_uv_update */
/* ========================================================================== */

/* Initialize the uv_update data structure */
void PIRO_BAND_LONG_NAME(init_uv_update)
(
    Int m,                          /* # rows in A                   */
    Int n,                          /* #columns in A                 */
    Int bl,                         /* lower bandwidth of A          */
    Int bu,                         /* upper bandwidth of A          */
    Int nr,                         /* # rows in the block           */
    PIRO_BAND_LONG_NAME(uv_update) *upd  /* data structure for the update */
)
{
    /* Assign initial bandwidth. This will not change as bl and bu changes */
    upd->bw = bl + bu ;

    upd->first = 1 ;

    if (nr > 1 && bl > 0 && bu > 1)
    {
        /* Unsymmetric matrix with nr > 1 */
        upd->t_blksize = bl - nr + 1 ;
        upd->blksize = bu - nr + 1 ;
    }
    else if (bu == 1)
    {
        /* lower Hessenberg matrix */
        upd->t_blksize = -1 ;
        upd->blksize = bl - nr + 2 ;
    }
    else if (bl == 0)
    {
        /* upper triangular unsymmetric matrix, or upper triangular part of
         * symmetric matrix */
        upd->t_blksize = -1 ;
        upd->blksize = bu - nr + 1 ;
    }
    else if (nr == 1)
    {
        /* Unsymmetric matrix with nr == 1 */
        upd->t_blksize = -1 ;
        upd->blksize = upd->bw ;
    }

    /* Assign the last row for the update of U and V */
    upd->vt_lrow_lower = MIN(upd->bw-1, n-1) ;
    if (nr == 1)
    {
        upd->vt_lrow_upper = MIN(upd->bw-1, n-1) ;
        upd->u_lrow_upper = MIN(upd->bw+bl-1, m-1) ;
    }
    else
    {
        upd->vt_lrow_upper = MIN(bu-1, n-1) ;
        upd->u_lrow_upper = MIN(upd->bw-1, m-1) ;
    }
    upd->u_lrow_lower = MIN(bl-1, m-1) ;

    /* Assign the first row for the update of U and V */
    upd->vt_srow_min = n-1 ;
    upd->u_srow_min = m-1 ;
}

/* ========================================================================== */
/* set_uv_update */
/* ========================================================================== */

/* Reset the uv_update data structure after reducing nr rows in the upper/lower
 * triangular part.
 * */
void PIRO_BAND_LONG_NAME(set_uv_update)
(
    Int m,                          /* # rows in A                   */
    Int n,                          /* #columns in A                 */
    Int bl,                         /* lower bandwidth of A          */
    Int bu,                         /* upper bandwidth of A          */
    Int nr,                         /* # rows in the block           */
    PIRO_BAND_LONG_NAME(uv_update) *upd  /* data structure for the uv update */
)
{
    Int swap ;
    Int incr ;

    if (upd->first == 1)
    {
        upd->first = 0 ;

        upd->vt_lrow_lower = MIN(upd->vt_lrow_lower+nr, n-1) ;
        upd->u_lrow_lower = MIN(upd->u_lrow_lower+nr, m-1) ;
        upd->vt_lrow_upper = MIN(upd->vt_lrow_upper+nr, n-1) ;
        upd->u_lrow_upper = MIN(upd->u_lrow_upper+nr, m-1) ;

        if (bl != 0 && nr > 1 && bu != 1)
        {
            /* Unsymmetric matrix with nr > 1 */
            SWAP(upd->t_blksize, upd->blksize, swap) ;
        }
    }
    else
    {
        incr = upd->blksize + nr - 1 ;
        upd->vt_lrow_upper = MIN(upd->vt_lrow_upper+incr , n-1) ;
        upd->u_lrow_upper = MIN(upd->u_lrow_upper+incr , m-1) ;

        if (bl != 0 && nr > 1 && bu != 1)
        {
            /* Unsymmetric matrix with nr > 1 */
            incr = upd->t_blksize + nr - 1 ;
        }

        upd->vt_lrow_lower = MIN(upd->vt_lrow_lower+incr, n-1) ;
        upd->u_lrow_lower = MIN(upd->u_lrow_lower+incr, m-1) ;
    }
}
