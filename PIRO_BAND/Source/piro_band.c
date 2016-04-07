/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band.c ========================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Reduces an m-by-n band matrix A to an upper bidiagonal matrix using blocked,
 * pipelined Givens rotations. The transformation is A = U * B * V'.
 *
 * The algorithm zeros out a block of entries nr*nc (nrl*ncl in lower band) and
 * chases the fill down the matrix. The pipelining of the rotations during the
 * chase will ensure that the fill during the chase will be no more than a
 * single scalar.
 *
 * The bidiagonal matrix will be returned as two vectors B1 and B2. The
 * rotations can also be accumulated as U and V.
 *
 * The rotations are accumulated and returned as V and not V'.

 TODO replace MALLOC and FREE with SuiteSparse_malloc and _free.
 TODO delete #include "piro_band_memory.h"
 * */

#include "piro_band_memory.h"
#include "piro_band_lapack_internal.h"
#include "piro_band_uv_update.h"
#include "piro_band.h"

#ifdef BENCHMARK
/* for internal benchmarking purposes only, not for normal usage */
double piro_band_flops = 0.0 ;
#endif

/* ========================================================================== */
/* === piro_band_reduce ===================================================== */
/* ========================================================================== */

int PIRO_BAND(reduce)
(
    Int blks[],
    /* An array of size four for the block size in upper and lower bands.
     * Specifically :
     * blk[0], blk[1] for #columns, #rows in the block for upper band
     * blk[2], blk[3] for #rows, #columns in the block for lower band
     * Block sizes cannot be negative. They can be zero only if the
     * corresponding bandwidth is zero.
     * blk[0] + blk[1] <= upper bandwidth
     * blk[2] + blk[3] <= lower bandwidth + 1 if lower bandwidth > 0
     * blk[2] + blk[3] = 0 otherwise.
     * */
    Int m,                  /* #rows in the original matrix.  m >= 0        */

    Int n,                  /* #columns in the original matrix. n >= 0      */

    Int nrc,                /* #rows in C
                             * nrc >= 0, and nrc = 0 if C == NULL           */

    Int bl,                 /* lower bandwidth,
                             * bl >= 0 and bl <= m-1 if m > 0.              */

    Int bu,                 /* upper bandwidth,
                             * bu >= 1 and bu <= n-1 if n > 0.              */

    Entry *A,               /* Band Matrix stored in packed band format,
    * If the matrix is complex then the real and imaginary parts should be
    * stored next to each other like C99 complex data type. zomplex where we
    * store real and imaginary parts are stored in two separate arrays are
    * not supported.
    * */

    Int ldab,               /* leading dimension of A
                             * ldab >= (bl + bu + 1)
                             * If A is symmetric then ldab >= bu + 1        */

    Entry *B1,              /* Output : diagonal 0 of the bidiagonal matrix
                             * A vector of size min(m, n). Always real      */

    Entry *B2,              /* Output : diagonal 1 of the bidiagonal matrix
                             * A vector of size min(m, n)-1. Always real.   */

    Entry *U,               /* Output : accumulated left rotations, mxm
                             * Orthogonal matrix. The rotations are always
                             * accumulated, and U always initialized to
                             * identity if U != NULL                        */

    Int ldu,                /* leading dimension of U
                             * ldu >= 0, ldu >= max(1, m) if U != NULL      */

    Entry *V,               /* Output :  accumulated right rotations, nxn
                             * orthogonal matrix. The rotations are always
                             * accumulated and V always initialized
                             * to identity if V != NULL.                    */

    Int ldv,                /* leading dimension of V
                             * ldv >= 0, ldv >= max(1, n) if V != NULL      */

    Entry *C,               /* Output : for left rotations  nrc x m matrix
                             * which accumulates the left rotations. C is
                             * never initialized.
                             * */

    Int ldc,                /* leading dimension of C, nrc >= 0, ldc >= nrc */

    Entry *dws,             /* workspace, to store the blocked Givens
    * rotations, of size 2*MAX (blk[0]*blk[1], blk[2]*blk[3]). If the input is
    * complex then we need twice this space to store the complex rotations.
    * In the complex case, We store both cosine and sine of the rotations as
    * complex scalars even though the cosine is always real for simpler
    * indexing.                                                             */

    Int sym                 /* 0: A is unsymmetric, nonzero: A is symmetric */
)
{
    int info ;
    Int iblk[4] ;
    Int msize ;
    Int esfmn2 ;
    Int alloc_work ;
    Entry safemin, eps ;
    PIRO_BAND_Common Common ;

    /* ----------------------------------------------------------------------*/
    /* Check the inputs  */
    /* ----------------------------------------------------------------------*/

    if (C == NULL)
    {
        /* ignore ldc and nrc */
        ldc = 0 ;
        nrc = 0 ;
    }
    if (U == NULL)
    {
        /* ignore ldu */
        ldu = 0 ;
    }
    if (V == NULL)
    {
        /* ignore ldv */
        ldv = 0 ;
    }
    if (sym != 0)
    {
        /* Treat i/p as symmetric */
        sym = 1;
    }

    info = PIRO_BAND(check) (m, n, nrc, bl, bu, A, ldab, B1, B2, U, ldu,
        V, ldv, C, ldc, sym, 1) ;
    if (info != PIRO_BAND_OK) return (info) ;
    if (m == 0 || n == 0) return (PIRO_BAND_OK) ;

    if (blks == NULL)
    {
        /* Ignore the workspace provided by the user */
        dws = NULL ;

        /* Get the recommended block size */
        PIRO_BAND_LONG_NAME(get_blocksize)(m, n, bl, bu,
            (C != NULL || U != NULL || V != NULL), iblk) ;
    }
    else
    {
        iblk[0] = blks[0] ;
        iblk[1] = blks[1] ;
        iblk[2] = blks[2] ;
        iblk[3] = blks[3] ;
    }

    /* Make sure block size is valid. Block size can't be negative or more than
     * the bandwidth or zero when there are entries to reduce
     * */
    if ( iblk[0] < 0 || iblk[1] < 0 || iblk[2] < 0 || iblk[3] < 0 ||
         (iblk[0] + iblk[1] > bu) || (bl > 0 && iblk[2] + iblk[3] > bl+1) ||
         (bl == 0 && iblk[2] + iblk[3] > 0) ||
         ((iblk[0] == 0 || iblk[1] == 0) && bu > 1) ||
         ((iblk[2] == 0 || iblk[3] == 0) && bl > 0) )
    {
        return (PIRO_BAND_BLKSIZE_INVALID) ;
    }

    /* ----------------------------------------------------------------------*/
    /* allocate workspace */
    /* ----------------------------------------------------------------------*/

    alloc_work = 0 ;
    if (dws == NULL && (bl > 0 || bu > 1))
    {
        /* Allocate workspace for saving the rotations */
        msize = CSIZE(2 * MAX (iblk[0] * iblk[1], iblk[2]*iblk[3])
                            * sizeof(Entry)) ;
        dws = (Entry *) MALLOC(msize) ;
        alloc_work = 1 ;
        if (dws == NULL)
        {
            return (PIRO_BAND_OUT_OF_MEMORY) ;
        }
    }

    /* ----------------------------------------------------------------------*/
    /* Initialize U and V */
    /* ----------------------------------------------------------------------*/

    if (U != NULL) PIRO_BAND(set_to_eye)(ldu, m, U) ;

    if (V != NULL) PIRO_BAND(set_to_eye)(ldv, n, V) ;

    /* ----------------------------------------------------------------------*/
    /* compute floating-point parameters for this computer */
    /* ----------------------------------------------------------------------*/

    /* Set the safe minimum and safe maximum */
    safemin = PIRO_BAND_Entry_MIN ;
    eps = PIRO_BAND_EPS/2 ;
    /* logf and powf are C99 requirements. Use log and pow always */
    esfmn2  = log(safemin/eps) / log(2) / 2 ;
    Common.safemin2 = pow(2, esfmn2) ;
    Common.safemx2 = 1/Common.safemin2 ;

    /* ----------------------------------------------------------------------*/
    /* reduce to bidiagonal form */
    /* ----------------------------------------------------------------------*/

    /* Call equal bandwidth code if it can be safely called. */
    if ((iblk[1] == iblk[3]) || sym || bl == 0  || bu == 1)
    {
        PIRO_BAND(reduce_equalbw)(A, ldab, m, n, nrc, bl, bu, iblk[0], iblk[1],
                iblk[2], iblk[3], dws, B1, B2, U, ldu, V, ldv, C, ldc, sym,
                &Common) ;
    }
    else
    {
        PIRO_BAND(reduce_unequalbw)(A, ldab, m, n, nrc, bl, bu, iblk[0],
                iblk[1], iblk[2], iblk[3], dws, B1, B2, U, ldu, V, ldv, C,
                ldc, sym, &Common) ;
    }

    if (alloc_work && dws != NULL)
    {
        FREE(dws) ;
    }
    return (PIRO_BAND_OK) ;
}

/* ========================================================================== */
/* === piro_band_reduce_equalbw ============================================= */
/* ========================================================================== */

/* Handle the case when upper bandwidth and the lower bandwidth are the same.
 *
 * When reducing the bandwidth this method alternates between the upper and
 * lower bands for each set of nr rows/columns. The algorithm for reduction is
 *
 *  1. While upper bandwidth > 1 or lower bandwidth > 0
 *      a. for each set of nr rows/columns
 *          reduce the lower bandwidth to nr-1 reducing nc rows in each column i
 *          in one iteration.
 *          reduce the upper bandwidth to nr reducing nc columns in each row in
 *          one iteration.
 *      b. Set new lower and upper bandwidth. Adjust block sizes such that new
 *      block will reduce entire column(row) in lower(upper) bandwidths to the
 *      final output (1xnr or 1xnrl).
 *
 * Note : Step b in the algorithm sets new block sizes to be 1xnr, because the
 * number of rows(columns) left in lower(upper) bandwidth is equal to nr which
 * will be much smaller than the original bandwidth. So all it can be reduced
 * as one block. So the while loop in step 1 will only be taken at most 2 times.
 *
 * For eg. If A is a 10x10 matrix with both upper and lower bandwidth equals 3
 * and the block size for the upper and lower methods equal to their respective
 * widths then the first two steps will be as shown below.
 *
 *   Initial Matrix           Step 1                      Step 2
 *   x x x x                   x x x x                    x x 0 0
 *   x x x x x                 0 x x x x                    x x x x
 *   x x x x x x               0 x x x x x                  x x x x x
 *   x x x x x x x             0 x x x x x x                x x x x x x
 *     x x x x x x x             x x x x x x x              x x x x x x x
 *       x x x x x x x             x x x x x x x              x x x x x x x
 *         x x x x x x x             x x x x x x x              x x x x x x x
 *           x x x x x x               x x x x x x                x x x x x x
 *             x x x x x                 x x x x x                  x x x x x
 *               x x x x                   x x x x                    x x x x
 *
 * */
void PIRO_BAND(reduce_equalbw)
(
    Entry *A,               /* Band Matrix                             */
    Int ldab,               /* leading dimension of A                  */
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    Int nc,                 /* #columns in the upper block             */
    Int nr,                 /* #rows in the upper block                */
    Int ncl,                /* #columns in the lower block             */
    Int nrl,                /* #rows in the lower block                */
    Entry *dws,             /* workspace of size 2*MAX(nc*nr, ncl*nrl) */
    Entry *B1,              /* o/p diagonal 0                          */
    Entry *B2,              /* o/p diagonal 1                          */
    Entry *U,               /* o/p accumulated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumulated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym,                /* Is the matrix symmetric ?               */
    PIRO_BAND_Common *Common /* common data structure                   */
)
{
    Int k ;
    Int i ;
    Int minmn ;
    Int obu, orig_bl ;
    Int m_nr ;
    Int adj ;
    Int swap ;
    Int want_uv ;                   /* 0/1 for need U or V or not       */
    /* temp double variables to find and apply Givens rotations */
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)] ;
    Entry d[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)], da1[CSIZE(1)], db1[CSIZE(1)] ;
    Entry da2[CSIZE(1)], db2[CSIZE(1)], dtemp2[CSIZE(1)] ;
    Entry con_s[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)], conjs1[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif
    PIRO_BAND_LONG_NAME(uv_update) update ;
    PIRO_BAND_LONG_NAME(uv_update) *upd ;
    upd = &update ;

    obu = bu ;
    orig_bl = bl ;
    minmn = MIN (m, n) ;
    m_nr = MAX (nr, nrl) ;

    /* ----------------------------------------------------------------------*/
    /* Initialization for the update of U and V */
    /* ----------------------------------------------------------------------*/

    want_uv = 0 ;
    if (U != NULL || V != NULL)
    {
        want_uv = 1 ;
        PIRO_BAND_LONG_NAME(init_uv_update)(m, n, bl, bu, m_nr, upd) ;
    }

    /* ----------------------------------------------------------------------*/
    /* Reduce to bidiagonal in two steps */
    /* ----------------------------------------------------------------------*/

    PRINT(("Reduce Equal bandwidths\n")) ;
    while (bl > 0 || bu > 1) /* This loop runs at most 2 times */
    {
        for (k = 0 ; k <= minmn ; k += m_nr)
        {
            upd->k = k ;
            ASSERT (m_nr > 0) ;
            /* Zero blocks of entries from columns k..k+m_nr-1  */
            if (!sym && bl > 0) /* bl > 0 not really required */
            {
                PIRO_BAND(reduce_blk_lower_band)(A, ldab, m, n, nrc, bl, bu,
                ncl, nrl, dws, k, obu, U, ldu, V, ldv, C, ldc, upd, Common) ;
            }

            if (want_uv)
            {
                /* Reset the update data structure */
                if (bl != 0 && m_nr > 1 && bu != 1)
                {
                    SWAP(upd->t_blksize, upd->blksize, swap) ;
                }
            }

            /* Zero blocks of entries from rows k..k+m_nr-1  */
            if (bu > 1)
            {
                PIRO_BAND(reduce_blk_upper_band)(A, ldab, m, n, nrc, bl, bu,
                    nc, nr, dws, k, obu, U, ldu, V, ldv, C, ldc, sym, upd,
                    Common) ;
            }

            if (want_uv)
            {
                /* Reset the update data structure */
                PIRO_BAND_LONG_NAME(set_uv_update)(m, n, bl, bu, m_nr, upd) ;
            }
        }

        adj = ((nrl == 1 || bu == 1) ? 0 : 1) ;
        if (want_uv)
        {
            /* U and V are almost full. */
            upd->vt_srow_min = m_nr ;
            upd->u_srow_min = m_nr-1+adj ;
            upd->blksize = upd->bw ;
            upd->vt_lrow_upper = 0 ;
            upd->vt_lrow_lower = 0 ;
            upd->u_lrow_lower = -bu ;
            upd->u_lrow_upper = bl ;
        }

        /* Adjust the lower bandwidth to the reduced size */
        /*bl = (bl == ncl) ? 0 : MIN (bl, nrl-1) ;*/
        bl = (nrl == 1) ? 0 : MIN (bl, nrl-1+adj) ;

        if (nr != 0) bu = MIN (bu, nr) ;
        PRINT(("bl="ID", bu="ID"\n", bl, bu)) ;

        nc = nr - 1 ;
        ncl = nrl - 1+adj ;
        /* set number of rows in blocks to 1, as band size is equal to original
         * block size now. For smaller bands we can reduce entire rows/columns
         * in one iteration effectively.
         * */
        nr = 1 ;
        nrl = 1 ;
        m_nr = 1 ;
    }

    /* ----------------------------------------------------------------------*/
    /* Reduce the entry in A(m-1, m) and chase the fill. */
    /* ----------------------------------------------------------------------*/

    if (m < n)
    {
        /* Chase the entry in A(m-1, m) */
        ASSIGN_TO_SCALAR(db, A, INDEX(m-1, m)) ;
        for (k = m-1 ; k >= 0 ; k--)
        {
            ASSIGN_TO_SCALAR(da, A, INDEX(k, k)) ;
            GIVENS(da, db, d, c, s, Common) ;
            ASSIGN_TO_MATRIX(d, A, INDEX(k, k)) ;
            if (k > 0)
            {
                /* Using da and d as temp */
                ASSIGN_TO_SCALAR(da, A, INDEX(k-1, k)) ;
                CONJ(con_s, s) ;
                NEG(con_s) ;
                MULT(db, con_s, da) ;
                MULT(d, c, da) ;
                ASSIGN_TO_MATRIX(d, A, INDEX(k-1, k)) ;
            }
            if (V != NULL)
            {
                APPLY_GIVENS_TO_COL_M(V, ldv, k, m, 0, n-1, c, s, i, da1, db1,
                                    dtemp, da2, db2, dtemp2, n) ;
            }
        }
    }

    /* ----------------------------------------------------------------------*/
    /* Convert the diagonal and the off diagonal elements to real */
    /* ----------------------------------------------------------------------*/

#ifdef PIROBAND_COMPLEX
    PIRO_BAND(complex_to_real)(A, ldab, m, n, nrc, orig_bl, obu, U, ldu, V, ldv,
                    C, ldc, sym ) ;
#endif

    /* ----------------------------------------------------------------------*/
    /* Copy the diagonals 0 and 1 to the output */
    /* ----------------------------------------------------------------------*/

    ASSIGN_MATRIX_TO_MATRIX(B1, CSIZE(0), A, CSIZE(obu)) ;
    i = obu ;
    PRINT(("B1[0] = "DID" \n", A[obu])) ;
    for (k = 1 ; k < minmn ; k++)
    {
        i += ldab ;
        PRINT(("B1["ID"] = "DID", B2["ID"] = "DID"\n", k, A[i], k-1, A[i-1])) ;
        /*Copy only the real part */
        B1[k] = A[CSIZE(i)] ;
        B2[k-1] = A[CSIZE(i-1)] ;
    }
}

/* ========================================================================== */
/* === piro_band_reduce_unequalbw =========================================== */
/* ========================================================================== */

/* Handles the case when upper bandwidth and the lower bandwidth aren't the same
 *
 * When reducing the bandwidth reduce the lower band fully before the super
 * diagonals. The Algorithn for the reduction is
 *
 *  1. While upper bandwidth > 1 or lower bandwidth > 0
 *      a. for each set of nr columns
 *          reduce the lower bandwidth to nr-1 reducing nc rows in each column
 *          in one iteration.
 *      b. Adjust the lower bandwidth.
 *      c. for each set of nr rows
 *          reduce the upper bandwidth to nr reducing nc columns in each row in
 *          one iteration.
 *      d. Set new upper bandwidth. Adjust block sizes such that new block will
 *       reduce entire column(row) in lower(upper) bandwidths to the final
 *       output (1xnr or 1xnrl).
 *
 * Note : Step d in the algorithm sets new block sizes to be 1xnr, because the
 * number of rows(columns) left in lower(upper) bandwidth is equal to nr which
 * will be much smaller than the original bandwidth. So all it can be reduced
 * as one block. So the while loop in step 1 will only be taken at most 2 times.
 *
 * For eg. If A is a 10x10 matrix with both upper and lower bandwidth equals 3
 * and 4 respectively and the block size for the upper and lower methods equal
 * to their respective widths then the first two steps will be as below.
 *
 *   Initial Matrix           Step 1                      Step 2
 *   x x x x                   x x x x                    x x x x
 *   x x x x x                 0 x x x x                    x x x x
 *   x x x x x x               0 x x x x x                  0 x x x x
 *   x x x x x x x             0 x x x x x x                0 x x x x x
 *   x x x x x x x x           0 x x x x x x x              0 x x x x x x
 *     x x x x x x x x           x x x x x x x x            0 x x x x x x x
 *       x x x x x x x x           x x x x x x x x            x x x x x x x x
 *         x x x x x x x             x x x x x x x              x x x x x x x
 *           x x x x x x               x x x x x x                x x x x x x
 *             x x x x x                 x x x x x                  x x x x x
 *
 * If the upper and lower bandwidths are not balanced enough choose different
 * block sizes for them.
 * */
void PIRO_BAND(reduce_unequalbw)
(
    Entry *A,               /* Band Matrix                             */
    Int ldab,               /* leading dimension of A                  */
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    Int nc,                 /* #columns in the upper block             */
    Int nr,                 /* #rows in the upper block                */
    Int ncl,                /* #columns in the lower block             */
    Int nrl,                /* #rows in the lower block                */
    Entry *dws,             /* workspace of size 2*MAX(nc*nr, ncl*nrl) */
    Entry *B1,              /* o/p diagonal 0                          */
    Entry *B2,              /* o/p diagonal 1                          */
    Entry *U,               /* o/p accumulated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumulated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym,                /* Is the matrix symmetric ?               */
    PIRO_BAND_Common *Common /* common data structure                   */
)
{
    Int k ;
    Int i ;
    Int minmn ;
    Int obu, orig_bl ;
    Int adj ;
    /* temp double variables to find and apply Givens rotations */
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)] ;
    Entry d[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)], da1[CSIZE(1)], db1[CSIZE(1)] ;
    Entry da2[CSIZE(1)], db2[CSIZE(1)], dtemp2[CSIZE(1)] ;
    Entry con_s[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)], conjs1[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif
    PIRO_BAND_LONG_NAME(uv_update) update ;
    PIRO_BAND_LONG_NAME(uv_update) *upd ;
    upd = &update ;

    orig_bl = bl ;
    obu = bu ;

    /* ----------------------------------------------------------------------*/
    /* Initialization for the update of U and V */
    /* ----------------------------------------------------------------------*/

    if (U != NULL || V != NULL)
    {
        /* Update the entire matrix. This will be slower */
        upd->bw = bl + bu ;
        upd->first = 1 ;
        upd->blksize = 1 ; /* No need for this, but just in case */
        /* Assign the last row for the update of U and V */
        upd->vt_lrow_lower = n-1 ;
        upd->vt_lrow_upper = n-1 ;
        upd->u_lrow_upper =  m-1 ;
        upd->u_lrow_lower = m-1 ;

        /* Assign the first row for the update of U and V */
        upd->vt_srow_min = n-1 ;
        upd->u_srow_min = m-1 ;
    }

    /* ----------------------------------------------------------------------*/
    /* Reduce the matrix to bidiagonal form in two steps */
    /* ----------------------------------------------------------------------*/

    minmn = MIN (m, n) ;
    PRINT(("Reduce different bandwidths\n")) ;
    while (bl > 0 || bu > 1) /* This loop runs at most 2 times */
    {
        for (k = 0 ; k <= minmn && !sym && bl > 0 ; k += nrl)
        {
            ASSERT (nrl > 0) ;
            upd->k = k ;
            /* Zero blocks of entries from columns k..k+nr-1  */
            PIRO_BAND(reduce_blk_lower_band)(A, ldab, m, n, nrc, bl, bu, ncl,
                    nrl, dws, k, obu, U, ldu, V, ldv, C, ldc, upd, Common) ;
        }

        /* Assign the first row for the update of U and V */
        upd->vt_srow_min = 0 ;
        upd->u_srow_min = 0 ;

        adj = ((nrl == 1 || bu == 1) ? 0 : 1) ;
        /* Adjust the lower bandwidth to the reduced size */

        if (nrl == 1)
        {
            bl = 0 ;
        }
        else
        {
            /* This (nrl != 1) case is not tested by Tcov tests.
               It is conjectured to be unreachable in any case.  We have it
               here just to be safe, in case our conjecture is incorrect. */
            bl = MIN (bl, nrl-1+adj) ;
        }

        for (k = 0 ; k <= minmn && bu > 1 ; k += nr)
        {
            ASSERT (nr > 0) ;
            upd->k = k ;
            /* Zero blocks of entries from rows k..k+nr-1  */
            PIRO_BAND(reduce_blk_upper_band)(A, ldab, m, n, nrc, bl, bu, nc, nr,
                    dws, k, obu, U, ldu, V, ldv, C, ldc, sym, upd, Common) ;
        }

        bu = MIN (bu, nr) ;
        PRINT(("bl="ID", bu="ID"\n", bl, bu)) ;

        nc = nr - 1 ;
        ncl = nrl - 1+adj ;
        /* set number of rows in block to 1, as band size is equal to original
         * block size now. For smaller bands we can reduce entire rows/columns
         * in one iteration effectively.
         * */
        nr = 1 ;
        nrl = 1 ;
    }

    /* ----------------------------------------------------------------------*/
    /* Reduce the entry in A(m-1, m) and chase the fill. */
    /* ----------------------------------------------------------------------*/

    if (m < n)
    {
        ASSIGN_TO_SCALAR(db, A, INDEX(m-1, m)) ;
        for (k = m-1 ; k >= 0 ; k--)
        {
            ASSIGN_TO_SCALAR(da, A, INDEX(k, k))
            GIVENS(da, db, d, c, s, Common) ;
            ASSIGN_TO_MATRIX(d, A, INDEX(k, k)) ;
            if (k > 0)
            {
                /* Using da and d as temp */
                ASSIGN_TO_SCALAR(da, A, INDEX(k-1, k)) ;
                CONJ(con_s, s) ;
                NEG(con_s) ;
                MULT(db, con_s, da) ;
                MULT(d, c, da) ;
                ASSIGN_TO_MATRIX(d, A, INDEX(k-1, k)) ;
            }
            if (V != NULL)
            {
                APPLY_GIVENS_TO_COL_M(V, ldv, k, m, 0, n-1, c, s, i, da1, db1,
                                dtemp, da2, db2, dtemp2, n) ;
            }
        }
    }

    /* ----------------------------------------------------------------------*/
    /* Convert the diagonal and the off diagonal elements to real */
    /* ----------------------------------------------------------------------*/

#ifdef PIROBAND_COMPLEX
    PIRO_BAND(complex_to_real)(A, ldab, m, n, nrc, orig_bl, obu, U, ldu, V, ldv,
                    C, ldc, sym ) ;
#endif

    /* ----------------------------------------------------------------------*/
    /* Copy the diagonals 0 and 1 to the output */
    /* ----------------------------------------------------------------------*/

    ASSIGN_MATRIX_TO_MATRIX(B1, CSIZE(0), A, CSIZE(obu)) ;
    i = obu ;
    PRINT(("B1[0] = "DID"\n ", A[obu])) ;
    for (k = 1 ; k < minmn ; k++)
    {
        i += ldab ;
        PRINT(("B1["ID"] = "DID", B2["ID"] = "DID"\n", k, A[i], k-1, A[i-1])) ;
        /*Copy only the real part */
        B1[k] = A[CSIZE(i)] ;
        B2[k-1] = A[CSIZE(i-1)] ;
    }
}

/* ========================================================================== */
/* === reduce_blk_lower_band ================================================ */
/* ========================================================================== */

/*
 * Zeros ncl*nrl entries as one block from the kth column.
 * Divides the matrix into blocks to zero the entries in seed block.
 * Chases the fill into other blocks till the end of the matrix.
 *
 *  There are five types of blocks
 * 1. S (seed) block, which covers the ncl*nrl entries to zero and other entries
 * to be touched in those columns. (marked z and S respectively)
 * 2. R block - row rotation block has no fill in it and only row rotations.
 * 3. F block - Fill block, uses both row and rotations to chase the fill. The
 * fill is no more than the tradition schwarz's method.
 * 4. C block - column rotation block has no fill in it and only column
 * rotations.
 * 5. D block - (lower) diagonal block. uses both row and column rotations.
 *
 * To zero the entries marked z, the blocks are split as below. The entries in
 * all the blocks are marked with the letters of their block names.
 * rseed is the pivotal row around which the nr*nc to be zeroed are defined.
 * The start/end row of a block is named srow and erow preceded by the block
 * name. For eg : dblk-srow is d-block start row.
 * The number that follows name (as d-blk srow1) is the iteration number.
 *
 *                               blk-ecol1      blk-ecol2
 *                     blk-scol1  |    blk-scol2 |
 *                             |  |        |     |
 *                             X  X  X  X
 *                             S  S  r  r  F  #                --> dblk-srow1
 *                             z  S  r  r  F  F  #
 *                             z  z  r  r  F  F  F  #          --> rseed1
 *                                z  r  r  F  F  F  F          --> dblk-erow1
 *                                   X  X  c  c  c  c  X
 *                                      X  c  c  c  c  X  X
 *                                         D  D  D  D  r  r    --> dblk-srow2
 *                                         #  D  D  D  r  r
 *                                            #  D  D  r  r    --> dblk-erow2
 *
 *   X - Original entries of matrix untouched in this iteration.
 *   z - Entries to be made zero.
 *   # - Fillin positions.(Not all at once, Fill in is never more than 1 scalar)
 *   S , r, F, c, D - Entries of S-block, R-block, F-block, C-Block and D-block
 *   respectively.
 *
 *   The Algorithm for the reduction :
 *   1. for each set of nc rows (in current nr columns starting at k)
 *      a. Divide the matrix into blocks
 *      b. Find row rotations to zero ncxnr entries in the S-block and apply
 *         them in in S-block.
 *      c. While there are row rotations
 *          c1. Apply row rotations to R-block.
 *          c2. Apply row rotations to F-block, Find column rotations to chase
 *              fill.
 *          c3. Adjust block start/end rows/columns.
 *          c3. Apply column rotations to C-block.
 *          c4. Apply column rotations to D-block, Find row rotations to chase
 *              fill.
 *
 * */
void PIRO_BAND(reduce_blk_lower_band)
(
    Entry *A,               /* Band Matrix                             */
    Int ldab,               /* leading dimension of A                  */
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    Int ncl,                /* #columns in the lower block             */
    Int nrl,                /* #rows in the lower block                */
    Entry *givens,          /* workspace for the rotations             */
    Int k,                  /* current column                          */
    Int obu,                /* orig. lower bandwidth of A, hidden usage :A() */
    Entry *U,               /* o/p accumulated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumulated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    PIRO_BAND_LONG_NAME(uv_update) *upd,  /* Structure for U and VT update */
    PIRO_BAND_Common *Common              /* common data structure         */
)
{
    Int iter ;                      /* iteration count           */
    Int cw ;                        /* current width             */
    Int nrow, ncol ;                /* # of rows/columns in an iteration */
    Int bscol, becol ;              /* block start/end column    */
    Int dsrow, derow ;              /* d-block start/end row     */
    Int rseed, cseed ;              /* row/column seeds          */
    Int csrow ;                     /* c-blk start row           */
    Int cgcnt ;                     /* column givens count       */
    Int rgcnt ;                     /* row givens count          */
    Int bw = 0 ;                    /* original bandwidth of A */
    Int vt_lrow = 0 ;               /* temp variable for the lastrow VT */
    Int u_lrow = 0 ;                /* temp variable for the lastrow of U */
    Int incr = 0 ;                  /* increment to the lastrow based on wave */
    Int want_uv ;                   /* 0/1 for need U or V or not       */
    Int adj ;

    iter = 0 ;
    cw = MIN (k+bl, m-1)-k  ;
    want_uv = 0 ;
    if (U != NULL || V != NULL)
    {
        want_uv = 1 ;
        bw = upd->bw ;
        incr = (upd->first == 0 ? upd->blksize : 1) ;
    }
    adj = ((nrl == 1 || bu == 1) ? 0 : 1) ;

    /* ----------------------------------------------------------------------*/
    /* Reduce one block of entries and changse the fill in lower band */
    /* ----------------------------------------------------------------------*/

    for ( ; cw >= nrl+adj ; cw -= ncl)
    {
        iter++ ;
        nrow = MIN (ncl, cw-(nrl-1+adj)) ;
        ncol = (k+nrl-1 > n-1) ? (n - k) : nrl ;

        /* create the blocks */
        bscol = k ;
        becol = MIN (k+ncol-1, n-1) ;
        dsrow = MAX (MIN (k+bl, m-1)-(ncl*iter), k+nrl-1+adj) ;
        derow = MIN (dsrow+nrow+ncol-1, m-1) ;
        rseed = MIN (dsrow+nrow, m-1) ;

        if (want_uv)
        {
            upd->c_lrow = MIN (upd->u_lrow_lower+incr, m-1) ;
        }
        /* Find the rotations from seed block */
        rgcnt = PIRO_BAND(sblk_lower_band)(m, nrc, rseed, dsrow+1, k, ncol,
            rseed, obu, ldab, A, givens, U, ldu, C, ldc, upd, Common) ;

        if (want_uv)
        {
            vt_lrow = upd->vt_lrow_lower ;
            u_lrow = MIN (upd->u_lrow_lower+bw, m-1) ;
        }
        while (1)
        {
            if (rgcnt == 0 || becol == n-1)
            {
                ASSERT (rgcnt == 0 || (ENDCOL(dsrow) == n-1)) ;
                break ;
            }

            /* Apply row rotations to R Block */
            PIRO_BAND(apply_rblk)(m, becol+1, ENDCOL(dsrow)-1, rgcnt, rseed,
                    dsrow, derow, obu, ldab, A, givens) ;

            cseed = ENDCOL(rseed) ;
            bscol = ENDCOL(dsrow) ;
            becol = ENDCOL(derow) ;

            PRINT(("bscol="ID", becol="ID"\n", bscol, becol)) ;
            if (bscol == becol)
            {
                /* This is just an optimization. F block will also handle it. */
                ASSERT (cseed == becol && cseed == n-1) ;
                PIRO_BAND(apply_rblk)(m, bscol, bscol, rgcnt, rseed, dsrow,
                    derow, obu, ldab, A, givens) ;
                break ;
            }

            /* Apply row rotations to F block, Get new column rotations */
            if (want_uv)
            {
                upd->c_lrow = MIN (vt_lrow+incr, n-1) ;
            }
            cgcnt = PIRO_BAND(apply_fblk)(m, n, rseed, dsrow+1, rgcnt, derow,
                        bscol, bu, obu, ldab, A, givens, V, ldv, upd, Common) ;

            if (derow == m-1 || cgcnt == 0)
            {
                /* No c blk or column rotations */
                break ;
            }

            rseed = MIN (cseed+bl, m-1) ;
            csrow = MIN (derow+1, m-1) ;
            dsrow = MIN (bscol+bl, m-1) ;
            derow = MIN (becol+bl, m-1) ;

            /* Apply column rotations to c blk */
            PIRO_BAND(apply_cblk)(n, cseed, bscol+1, cgcnt, csrow, dsrow,
                becol, obu, ldab, A, givens) ;

            /* Apply column rotations to D block, get new row rotations */
            if (want_uv)
            {
                upd->c_lrow = MIN (u_lrow+incr, m-1) ;
            }
            rgcnt = PIRO_BAND(apply_dblk)(m, n, nrc, cseed, bscol+1, cgcnt,
                dsrow, becol, bl, obu, ldab, A, givens, U, ldu, C, ldc,
                0, upd, Common) ;

            if (want_uv)
            {
                vt_lrow = MIN (vt_lrow+bw, n-1) ;
                u_lrow = MIN (u_lrow+bw, m-1) ;
            }
        }

    }
}

/* ========================================================================== */
/* === reduce_blk_upper_band ================================================ */
/* ========================================================================== */

/*
 * Zeros ncl*nrl entries as one block from the kth row.
 * Divides the matrix into blocks to zero the entries in seed block.
 * Chases the fill into other blocks till the end of the matrix.
 *
 * For details of the blocks, check the comments above reduce_blk_lower_band.
 * The picture below will illustrate the blocks to zero out 2x2 entries marked
 * z below.
 *
 * To zero the entries marked z, the blocks are split as below. The entries in
 * all the blocks are marked with the letters of their block names.
 * cseed is the pivotal column around which the nr*nc to be zeroed are defined.
 * The start/end row of a block is named srow and erow preceded by the block
 * name. For eg : dblk-srow is d-block start row.
 * The number that follows name (as d-blk srow1) is the iteration number.
 *
 *
 *                                  column-seed1
 *                           blk-scol1  |  blk-ecol1
 *                                |     |  |
 *                             X  S  z  z
 *                             X  S  S  z  z                  --> s-blk erow1
 *                             X  r  r  r  r  X               --> c-blk srow1
 *                             X  r  r  r  r  X  X
 *                                D  D  D  D  c  c  F  .      --> d-blk srow1
 *                                .  D  D  D  c  c  F  F  .
 *                                   .  D  D  c  c  F  F  F   --> row-seed1
 *                                      .  D  c  c  F  F  F   --> d-blk erow1
 *                                            X  X  r  r  r
 *                                               X  r  r  r
 *
 *
 *   X - Original entries of matrix untouched in this iteration.
 *   z - Entries to be made zero.
 *   # - Fillin positions.(Not all at once, Fill in is never more than 1 scalar)
 *   S , r, F, c, D - Entries of S-block, R-block, F-block, C-Block and D-block
 *   respectively.
 *
 *   The Algorithm for the reduction :
 *   1. for each set of nc columns (in current nr rows starting at k)
 *      a. Divide the matrix into blocks
 *      b. Find column rotations to zero ncxnr entries in the S-block and apply
 *         them in in S-block.
 *      c. While there are column rotations
 *          c1. Apply column rotations to C-block.
 *          c2. Apply column rotations to D-block, Find row rotations to chase
 *              fill.
 *          c3. Apply row rotations to R-block.
 *          c4. Apply row rotations to F-block, Find column rotations to chase
 *              fill.
 *          c5. Adjust block start/end rows/columns.
 *
 * The blocking technique for the upper band is the same as lower band except
 * that the order of the blocks in step 1c are different. C, D, R, F is the
 * sequence here. While lower band uses R, F, C, D.
 *
 * */
void PIRO_BAND(reduce_blk_upper_band)
(
    Entry *A,               /* Band Matrix                             */
    Int ldab,               /* leading dimension of A                  */
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    Int nc,                 /* #columns in the upper block             */
    Int nr,                 /* #rows in the upper block                */
    Entry *givens,          /* workspace for the rotations             */
    Int k,                  /* current row                             */
    Int obu,                /* orig. lower bandwidth of A, hidden usage :A() */
    Entry *U,               /* o/p accumulated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumulated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym,                /* Is the matrix symmetric ?               */
    PIRO_BAND_LONG_NAME(uv_update) *upd, /* Structure for U and VT update */
    PIRO_BAND_Common *Common             /* common data structure         */
)
{
    Int iter ;                      /* iteration count                */
    Int cw ;                        /* current width                  */
    Int ncol, nrow ;                /* # of cols in an iteration      */
    Int bscol, becol ;              /* block start/end column         */
    Int dsrow, derow ;              /* d-block start/end row          */
    Int rseed, cseed ;              /* row/column seeds               */
    Int csrow, serow ;              /* c-blk start row, s-blk end row */
    Int cgcnt ;                     /* column givens count            */
    Int rgcnt ;                     /* row givens count               */
    Int bw = 0 ;                    /* original bandwidth of A */
    Int u_lrow = 0 ;                /* temp variable for the lastrow U     */
    Int vt_lrow = 0 ;               /* temp variable for the lastrow of VT */
    Int incr = 0 ;                  /* increment to the lastrow based on wave */
    Int want_uv ;                   /* 0/1 for need U or V or nt        */

    iter = 0 ;
    cw = ENDCOL(k)-k ;
    want_uv = 0 ;
    if (U != NULL || V != NULL)
    {
        want_uv = 1 ;
        bw = upd->bw ;
        incr = (upd->first == 0 ? upd->blksize : 1) ;
    }

    /* ----------------------------------------------------------------------*/
    /* Reduce one block of entries and changse the fill in lower band */
    /* ----------------------------------------------------------------------*/
    for ( ; cw > nr ; cw -= nc)
    {
        iter++ ;
        ncol = MIN (nc, cw-nr) ;
        nrow = (k+nr-1 > m-1) ? (m - k) : nr ;

        /* Create blocks */
        bscol = MAX (ENDCOL(k)-(nc*iter), k+nr) ;
        becol = MIN (bscol+ncol+nrow-1, n-1) ;
        serow = MIN (k+nrow-1, m-1) ;
        dsrow = MIN (bscol+bl, m-1) ;
        derow = MIN (becol+bl, m-1) ;
        cseed = bscol + ncol ;
        csrow = MIN (serow+1, m-1) ;

        if (want_uv)
        {
            upd->c_lrow = MIN (upd->vt_lrow_upper+incr, n-1) ;
        }
        /* Find the rotations from seed block */
        cgcnt = PIRO_BAND(sblk_upper_band)(n, cseed, bscol+1, k, nrow, becol,
                            serow, cseed, obu, ldab, A, givens, V, ldv, upd,
                            Common) ;

        if (want_uv)
        {
            u_lrow = upd->u_lrow_upper ;
            vt_lrow = MIN (upd->vt_lrow_upper+bw, n-1) ;
        }
        while (1)
        {
            if (serow == m-1 || cgcnt == 0)
            {
                break ;
            }

            /* Apply column rotations to C blk */
            PIRO_BAND(apply_cblk)(n, cseed, bscol+1, cgcnt, csrow, dsrow, becol,
                        obu, ldab, A, givens) ;

            /* Apply column rotations to D block, get new row rotations */
            if (want_uv)
            {
                upd->c_lrow = MIN (u_lrow+incr, m-1) ;
            }
            rgcnt = PIRO_BAND(apply_dblk)(m, n, nrc, cseed, bscol+1, cgcnt,
                    dsrow, becol, bl, obu, ldab, A, givens, U, ldu, C, ldc, sym,
                    upd, Common);

            if (rgcnt == 0 || becol == n-1)
            {
                break ;
            }

            rseed = MIN (cseed+bl, m-1) ;

            /* Apply row rotations to R block */
            PIRO_BAND(apply_rblk)(m, becol+1, ENDCOL(dsrow)-1, rgcnt, rseed,
                    dsrow, derow, obu, ldab, A, givens) ;

            cseed = ENDCOL(MIN (cseed+bl, m-1)) ;
            bscol = ENDCOL(MIN (bscol+bl, m-1)) ;
            becol = ENDCOL(MIN (becol+bl, m-1)) ;

            if (bscol == becol)
            {
                /* This is just an optimization. F block will also handle it. */
                ASSERT (cseed == becol && cseed == n-1) ;
                PIRO_BAND(apply_rblk)(m, bscol, bscol, rgcnt, rseed, dsrow,
                    derow, obu, ldab, A, givens) ;
                break ;
            }

            /* Apply row rotations to F block, Get new column rotations */
            if (want_uv)
            {
                upd->c_lrow = MIN (vt_lrow+incr, n-1) ;
            }
            cgcnt = PIRO_BAND(apply_fblk)(m, n, rseed, dsrow+1, rgcnt, derow,
                        bscol, bu, obu, ldab, A, givens, V, ldv, upd, Common) ;

            if (derow == m-1)
            {
                /* No c blk */
                break ;
            }

            csrow = MIN (derow+1, m-1) ;
            dsrow = MIN (bscol+bl, m-1) ;
            derow = MIN (becol+bl, m-1) ;
            if (want_uv)
            {
                u_lrow = MIN (u_lrow+bw, m-1) ;
                vt_lrow = MIN (vt_lrow+bw, n-1) ;
            }

        } /* while 1 chase */
    }
}

/* ========================================================================== */
/* === sblk_upper_band ====================================================== */
/* ========================================================================== */

/* The S block is defined by rows k..k+nr-1 (serow) and columns  bscol..becol
 * Algorithm:
 * 1. For all rows in S-block
 *      a. If not the first row  apply column rotations from the previous rows.
 *      b. Find the column rotations to zero entries in current row (if there
 *      are any)
 *
 * Note : The block size is fixed in the number of rows. So we may have rows in
 * seed block, but may not be able find new rotations in them for the last few
 * rows of the matrix. We will simply apply rotations from the previous rows
 * then.
 *
 * The process of creating the rotations is like forming a wave. The first wave
 * from row k, passes row k+1, and a new wave of rotations are also generated
 * from k+1. Rotations from both k and k+1 will both be applied to row k+2
 * before forming its own rotations.
 *
 * Returns the number of new column rotations in the block.
 * */
Int PIRO_BAND(sblk_upper_band)
(
    Int n,              /* #columns in the matrix     */
    Int scol,           /* start column for the s-blk */
    Int ecol,           /* end column for the s-blk   */
    Int k,              /* current row                */
    Int nrow,           /* nr for the current block   */
    Int becol,          /* only for debug */
    Int serow,          /* only for debug */
    Int cseed,          /* column seed for the s-blk   */
    Int obu,            /* original lower bandwidth  of A, hidden usage :A() */
    Int ldab,           /* leading dim of the A, usage hidden in A() macro */
    Entry *A,           /* Band matrix */
    Entry *givens,      /* workspace for the rotations */
    Entry *V,           /* r.h.s of A to apply the rotations */
    Int ldv,            /* leading dimension of v */
    PIRO_BAND_LONG_NAME(uv_update) *upd,  /* Structure for U and VT update */
    PIRO_BAND_Common *Common              /* common data structure         */
)
{
    Int cgcnt ;
    Int row, col ;
    Int lcol ;
    Int tlcol ;
    Int ccol ;
    Int cnt ;
    Int tk ;
    Int i ;
    /* temproary variable for uv update */
    Int c_lrow = 0 ;
    Int c_srow = 0 ;
    Int s1 ;
    /* temp double variables to find and apply Givens rotations */
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)] ;
    Entry d[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)] ;
    Entry da1[CSIZE(1)], db1[CSIZE(1)], dtemp1[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)], conjs1[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif

    PRINT(("S block Upper Band\n")) ;
    cgcnt = 0 ;
    lcol = ecol ;
    /*Initialize for V update */
    if (V != NULL)
    {
        c_lrow = upd->c_lrow ;
        c_srow = upd->k + 1 ;
    }

    /* ----------------------------------------------------------------------*/
    /* Reduce one blocks of entries and save the rotations for the chase */
    /* ----------------------------------------------------------------------*/

    for (row = k ; row <= k+nrow-1 ; row++)
    {
        ASSERT (lcol <= becol && row <= serow) ;
        cnt = 0 ;
        tlcol = ecol ;
        col = cseed ;
        /* If not first row, Apply pending col rotations to current row */
        for (tk = k ; tk < row ; tk++, tlcol++)
        {
            PIRO_BAND_FLOPS((col-tlcol+1)*6) ;
            for (ccol = col ; ccol >= tlcol ; ccol--)
            {
                ASSIGN_CS_SCALARS(c, s, givens, cnt) ;
                cnt++ ;
                PRINT(("Apply crot="ID"\n", cnt)) ;
                PRINT(("row="ID" col="ID"-"ID"\n", row, ccol-1, ccol)) ;

                ASSIGN_TO_SCALAR(da, A, INDEX(row, ccol-1)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(row, ccol)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(row, ccol-1)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(row, ccol)) ;

            }
            col = MIN (col+1, n-1) ;
        }

        PIRO_BAND_FLOPS(6*(scol-lcol+1)) ;
        /* Find new col rotations to zero entries in current row */
        for (col = scol ; col >= lcol ; col--)
        {
            ASSIGN_TO_SCALAR(da, A, INDEX(row, col-1)) ;
            ASSIGN_TO_SCALAR(db, A, INDEX(row, col)) ;

            GIVENS(da, db, d, c, s, Common) ;
            SAVE_CS_SCALARS(givens, cgcnt, c, s) ;
            cgcnt++ ;
            PRINT(("find crot="ID"\n", cgcnt)) ;
            PRINT(("row="ID" col="ID"-"ID"\n", row, col-1, col)) ;

            ASSIGN_TO_MATRIX(d, A, INDEX(row, col-1)) ;
            ASSIGN_ZERO_TO_MATRIX(A, INDEX(row, col)) ;
            /* A[row, col] = 0.0 */

            if (V != NULL)
            {
                s1 = MIN (col-c_srow, upd->vt_srow_min) ;
                /* Need to return V, so apply the rotations to V immediately */
                APPLY_GIVENS_TO_COL_M(V, ldv, col-1, col, s1, c_lrow, c, s,
                                    i, da, db, dtemp, da1, db1, dtemp1, n) ;
            }
        }
        /* Break if end of the matrix, no more useful rotations available
         * in this seed block. */
        if (scol == lcol && scol == n-1)
            break ; /* Could remove this break and add lcol <= scol to
                    * first loop and the additions to loop incr */

        lcol++ ;
        scol = MIN (scol+1, n-1) ;

        if (V != NULL)
        {
            c_lrow = MIN (c_lrow+1, n-1) ;
            c_srow++ ;
        }
    }

    /* ----------------------------------------------------------------------*/
    /* If we did break at the end of the matrix, Apply all the column rotations
     * in the rest of the seed block. */
    /* ----------------------------------------------------------------------*/

    for (row = row+1 ; row <= k+nrow-1 ; row++)
    {
        cnt = 0 ;
        tlcol = ecol ;
        col = cseed ;
        /* Apply pending col rotations to current row */
        for (tk = k ; tk < row ; tk++, tlcol++)
        {
            PIRO_BAND_FLOPS(6*(col-tlcol+1)) ;
            for (ccol = col ; ccol >= tlcol ; ccol--)
            {
                ASSIGN_CS_SCALARS(c, s, givens, cnt) ;
                cnt++ ;
                PRINT(("Apply crot="ID"\n", cnt)) ;
                PRINT(("row="ID" col="ID"-"ID"\n", row, ccol-1, ccol)) ;

                ASSIGN_TO_SCALAR(da, A, INDEX(row, ccol-1)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(row, ccol)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(row, ccol-1)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(row, ccol)) ;
            }
            col = MIN (col+1, n-1) ;
        }
    }

    /* return the number of column rotations in this block */
    return cgcnt ;
}

/* ========================================================================== */
/* === sblk_lower_band ====================================================== */
/* ========================================================================== */

/* The S block is defined by columns k..k+nr-1 and rows  dsrow..derow
 * Algorithm:
 * 1. For all rows in S-block
 *      a. If not the first column  apply row rotations from the previous
 *      columns.
 *      b. Find the row rotations to zero entries in current column (if there
 *      are any)
 *
 * Note : The block size is fixed in the number of columns. So we may have
 * columns in seed block, but may not be able find new rotations in them for the
 * last few columns of the matrix. We will simply apply rotations from the
 * previous columns then.
 *
 * The process of creating the rotations is like forming a wave. The first wave
 * from column k, passes column k+1, and a new wave of rotations are also
 * generated from k+1. Rotations from both k and k+1 will both be applied to
 * column k+2 before forming its own rotations.
 *
 * Returns the number of new row rotations in the block.
 *
 * */
Int PIRO_BAND(sblk_lower_band)
(
    Int m,              /* #rows in the matrix      */
    Int nrc,            /* #columns in C            */
    Int srow,           /* start row of s-blk       */
    Int erow,           /* end row of s-blk         */
    Int k,              /* current column           */
    Int ncol,           /* #rows in the blk         */
    Int rseed,          /* row seed for the s-blk   */
    Int obu,            /* original lower b/w of A, usage hidden : A() macro */
    Int ldab,           /* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix                             */
    Entry *givens,      /* workspace for the rotations             */
    Entry *U,           /* o/p accumulated left rotations          */
    Int ldu,            /* leading dimension of u                  */
    Entry *C,           /* o/p for left rotations                  */
    Int ldc,            /* leading dimension of C                  */
    PIRO_BAND_LONG_NAME(uv_update) *upd, /* Structure for U and VT update */
    PIRO_BAND_Common *Common              /* common data structure         */
)
{
    Int rgcnt ;
    Int row, col ;
    Int lrow ;
    Int terow ;
    Int crow ;
    Int cnt ;
    Int tk ;
    Int i ;
    /* temporary variable for uv update */
    Int c_lrow = 0 ;
    Int c_srow = 0 ;
    Int s1 ;
    /* temp double variables to find and apply Givens rotations */
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)] ;
    Entry d[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)] ;
    Entry conj_s[2] ;
    Entry da1[CSIZE(1)], db1[CSIZE(1)], dtemp1[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)], conjs1[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif

    PRINT(("S block Lower Band\n")) ;
    rgcnt = 0 ;
    lrow = erow ;

    /* for U update */
    if (U != NULL)
    {
        c_lrow = upd->c_lrow ;
        c_srow = upd->k + 1 ;
    }

    /* ----------------------------------------------------------------------*/
    /* Reduce one blocks of entries and save the rotations for the chase */
    /* ----------------------------------------------------------------------*/

    for (col = k ; col <= k+ncol-1 ; col++)
    {
        cnt = 0 ;
        terow = lrow ;
        row = rseed ;
        /* Apply pending row rotations to current col */
        for (tk = k ; tk < col ; tk++, terow++)
        {
            PIRO_BAND_FLOPS(6*(row-terow+1)) ;
            for (crow = row ; crow >= terow ; crow--)
            {
                ASSIGN_CS_SCALARS(c, s, givens, cnt) ;
                cnt++ ;
                PRINT(("Apply rot="ID"\n", cnt)) ;
                PRINT(("row="ID"-"ID" col="ID"\n", crow-1, crow, col)) ;

                ASSIGN_TO_SCALAR(da, A, INDEX(crow-1, col)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(crow, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(crow-1, col)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(crow, col)) ;

            }
            row = MIN (row+1, m-1) ;
        }

        /* Find row rotations to zero entries in col */
        PIRO_BAND_FLOPS(6*(srow-erow+1)) ;
        for (row = srow ; row >= erow ; row--)
        {
            ASSIGN_TO_SCALAR(da, A, INDEX(row-1, col)) ;
            ASSIGN_TO_SCALAR(db, A, INDEX(row, col)) ;

            GIVENS(da, db, d, c, s, Common) ;

            SAVE_CS_SCALARS(givens, rgcnt, c, s) ;
            rgcnt++ ;
            PRINT(("find rot="ID"\n", rgcnt)) ;
            PRINT(("row="ID"-"ID" col="ID"\n", row-1, row, col)) ;

            ASSIGN_TO_MATRIX(d, A, INDEX(row-1, col)) ;
            ASSIGN_ZERO_TO_MATRIX(A, INDEX(row, col)) ;

            CONJ(conj_s, s) ;
            if (U != NULL)
            {
                s1 = MIN (row-c_srow, upd->u_srow_min) ;
                APPLY_GIVENS_TO_COL_M(U, ldu, row-1, row, s1, c_lrow, c, conj_s,
                                    i, da, db, dtemp, da1, db1, dtemp1, m) ;
            }
            if (C != NULL)
            {
                APPLY_GIVENS_TO_COL_M(C, ldc, row-1, row, 0, nrc-1, c, conj_s,
                                    i, da, db, dtemp, da1, db1, dtemp1, m) ;
            }
        }
        if (srow == erow && srow == m-1)
            break ;

        erow++ ;
        srow = MIN (srow+1, m-1) ;

        if (U != NULL)
        {
            c_lrow = MIN (c_lrow+1, m-1) ;
            c_srow++ ;
        }
    }

    /* ----------------------------------------------------------------------*/
    /* If we did break at the end of the matrix, Apply all the column rotations
     * in the rest of the seed block. */
    /* ----------------------------------------------------------------------*/

    PRINT(("col = "ID"\n", col)) ;
    for (col = col+1 ; col <= k+ncol-1 ; col++)
    {
        cnt = 0 ;
        terow = lrow ;
        row = rseed ;
        /* Apply pending row rotations to current col */
        for (tk = k ; tk < col ; tk++, terow++)
        {
            PIRO_BAND_FLOPS(6*(row-terow+1)) ;
            for (crow = row ; crow >= terow ; crow--)
            {
                ASSIGN_CS_SCALARS(c, s, givens, cnt) ;
                cnt++ ;
                PRINT(("Apply rot="ID"\n", cnt)) ;
                PRINT(("row="ID"-"ID" col="ID"\n", crow-1, crow, col)) ;

                ASSIGN_TO_SCALAR(da, A, INDEX(crow-1, col)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(crow, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(crow-1, col)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(crow, col)) ;
            }
            row = MIN (row+1, m-1) ;
        }
    }
    /* return the number of row rotations in this block */
    return rgcnt ;
}

/* ========================================================================== */
/* === apply_rblk =========================================================== */
/* ========================================================================== */

/* Applies rgcnt row rotations to the R-block. The block is defined by columns
 * scol..ecol and rows dsrow..derow. The rotations are applied centered around
 * rseed.
 *
 * If the R-block is a 6x2 submatrix and there are two waves each with 4 row
 * rotations be applied to the R-block, from a 2x4 block size, then applying
 * the rotations can be described pictorially as below. Rotations from the
 * first wave (G1..G4) will be applied to any column before rotations from the
 * second wave (G5..G8).
 *
 *                               R  R
 *                        G4
 *                               R  R
 *                   G8   G3
 *                               R  R
 *                   G7   G2
 *                               R  R
 *                   G6   G1
 *                               R  R     --> rseed, where first rotation starts
 *                   G5
 *                               R  R
 *  R - R-block entries.
 *  G1..G8 - Givens rotations to be applied to the R-block.
 *  Algorithm:
 *  1. for each column in R-block
 *      for each wave of rotations
 *          Apply all the rotations in current wave to appropriate rows in
 *          current column.
 *          Shift row indices for next wave.
 * */
void PIRO_BAND(apply_rblk)
(
    Int m,              /* #rows in the matrix              */
    Int scol,           /* starting column of R block       */
    Int ecol,           /* end column of R block            */
    Int rgcnt,          /* #rotations to be applied         */
    Int rseed,          /* row seed for the R block         */
    Int dsrow,          /* end row for the first wave       */
    Int derow,          /* used only for ASSERT             */
    Int obu,            /* original lower b/w of A, usage hidden : A() macro */
    Int ldab,           /* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix                             */
    Entry *givens       /* workspace for the rotations             */
)
{
    Int row, col ;
    Int srow, erow ;
    Int cnt ;
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif

    PRINT(("R block\n")) ;
    for (col = scol ; col <= ecol ; col++)
    {
        erow = dsrow + 1 ;
        srow = rseed ;
        cnt = 0 ;
        while (cnt < rgcnt)
        {
            ASSERT (srow >= erow && srow <= derow) ;
            PIRO_BAND_FLOPS(6*(srow-erow+1)) ;
            ASSIGN_TO_SCALAR(da, A, INDEX(srow, col)) ;
            for (row = srow ; row >= erow ; row--)
            {
                ASSIGN_CS_SCALARS(c, s, givens, cnt) ;
                cnt++ ;
                PRINT(("Apply rot="ID"\n", cnt)) ;
                PRINT(("row="ID"-"ID" col="ID"\n", row-1, row, col)) ;

                ASSIGN_TO_SCALAR(db, da, 0) ;
                ASSIGN_TO_SCALAR(da, A, INDEX(row-1, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(row, col)) ;
            }
            ASSIGN_TO_MATRIX(da, A, INDEX(erow-1, col)) ;

            erow++ ;
            srow = MIN (srow+1, m-1) ;
        }
    }
}

/* ========================================================================== */
/* === apply_cblk =========================================================== */
/* ========================================================================== */

/* Applies cgcnt column rotations to the R-block. The block is defined by
 * columns scol..ecol and rows srow..erow. The rotations are applied centered
 * around cseed.
 *
 * If the C-block is a 2x6 submatrix and there are two waves each with 4 row
 * rotations be applied to the C-block, from a 2x4 block size, then applying
 * the rotations can be described pictorially as below. Rotations from the
 * first wave (G1..G4) will be applied to any row before rotations from the
 * second wave (G5..G8).
 *
 *                            G8 G7 G6 G5
 *                         G4 G3 G2 G1
 *
 *                        C  C  C  C  C  C
 *                        C  C  C  C  C  C
 *  R - R-block entries.
 *  G1..G8 - Givens rotations to be applied to the C-block.
 *
 *  Though the pictorial representation describes it in terms of rows, the
 *  actual algorithm applies the rotation on columns.
 *
 *  Algorithm:
 *  1. for each wave of column rotation
 *      a. for each column in the appropriate set of columns
 *          Apply all the rotations in current wave to the rows in
 *          current column.
 *      b. Shift column indices for next wave.
 * */
void PIRO_BAND(apply_cblk)
(
    Int n,              /* #columns in the matrix           */
    Int scol,           /* start column for the first wave  */
    Int ecol,           /* end column for the first wave    */
    Int cgcnt,          /* #rotations to be applied         */
    Int srow,           /* start column for the C block     */
    Int erow,           /* end column for the C block       */
    Int becol,          /* used only for ASSERT             */
    Int obu,            /* original lower b/w of A, usage hidden : A() macro */
    Int ldab,           /* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix                      */
    Entry *givens       /* workspace for the rotations      */
)
{
    Int row, col ;
    Int cnt ;
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif

    PRINT(("C block\n")) ;
    if (srow >= erow)
        return ;

    for (cnt = 0 ; cnt < cgcnt ; )
    {
        ASSERT (scol >= ecol && scol <= becol) ;
        PIRO_BAND_FLOPS(6*(erow-srow)*(scol-ecol+1)) ;
        for (col = scol ; col >= ecol ; col--)
        {
            ASSIGN_CS_SCALARS(c, s, givens, cnt) ;
            cnt++ ;
            PRINT(("Apply crot="ID"\n", cnt)) ;
            PRINT(("row="ID"-"ID" col="ID"-"ID"\n", srow, erow-1, col-1, col));

            for (row = srow ; row < erow ; row++)
            {
                ASSIGN_TO_SCALAR(da, A, INDEX(row, col-1)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(row, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(row, col-1)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(row, col)) ;
            }
        }
        scol = MIN (scol+1, n-1) ;
        ecol++ ;
    }
}

/* ========================================================================== */
/* === apply_fblk =========================================================== */
/* ========================================================================== */

/* The F block is defined by rows dsrow..derow and columns bscol..ENDCOL(derow)
 * Apply row rotations, and remove any fill w/ column rotations and apply those.
 *
 * For every wave of row rotations that comes into the F block we will
 * need two passes on the F block.
 * Pass one will not generate any fill, but will apply all the row rotations to
 * the columns w/o any fill. This is a column operation even to apply the row
 * rotations. (Step 1a below)
 * Pass two will create the fill with row rotation and remove the fill using a
 * column rotation and will apply the column rotation immediately to the entire
 * column. (Step 1b below)
 *
 * Algorithm :
 *  1. For each wave of row rotations
 *      a. for each column in F-block (left-right)
 *          Apply all row rotations from current wave, except those that
 *          generate fillin in current column, to current column.
 *      b. for each column in F-block (right-left, that is the order for fill)
 *          Apply the pending rotation that will create fill.
 *          Remove fill with a column rotation.
 *          Apply column rotation to entire column.
 *      c. Adjust row indices for the next wave.
 *
 * Returns the number of new column rotations in the block.
 * */
Int PIRO_BAND(apply_fblk)
(
    Int m,              /* #rows in the matrix              */
    Int n,              /* #columns in the matrix           */
    Int srow,           /* start row for the first wave     */
    Int erow,           /* end row for the first wave       */
    Int rgcnt,          /* #rotations to be applied         */
    Int derow,          /* end row of the F block           */
    Int bscol,          /* start column of the F block      */
    Int bu,             /* upper bandwidth of the matrix    */
    Int obu,            /* original lower b/w of A, usage hidden : A() macro */
    Int ldab,           /* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix                      */
    Entry *givens,      /* workspace for the rotations      */
    Entry *V,           /* r.h.s of A to apply the rotations*/
    Int ldv,            /* leading dimension of v           */
    PIRO_BAND_LONG_NAME(uv_update) *upd, /* Structure for U and VT update */
    PIRO_BAND_Common *Common             /* common data structure         */
)
{
    Int cgcnt ;
    Int cnt, tcnt ;
    Int terow ;
    Int row, col ;
    Int crow ;
    Int i ;
    /* temproary variable for uv update */
    Int c_lrow = 0 ;
    Int c_srow = 0 ;
    Int s1 ;
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)], ws[CSIZE(1)] ;
    Entry d[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)] ;
    Entry da1[CSIZE(1)], db1[CSIZE(1)], dtemp1[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)], conjs1[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif

    PRINT(("F block\n")) ;
    cgcnt = 0 ;

    /* V update */
    if (V != NULL)
    {
        c_lrow = upd->c_lrow ;
        c_srow = upd->k + 1 ;
    }

    for (cnt = 0 ; cnt < rgcnt ; )
    {
        ASSERT (srow >= erow && srow <= derow) ;
        terow = erow ;

        /*-------------------- Phase 1 ----------------- */
        /*Apply the row rotation in current wave without generating any fill */
        for (col = bscol ; col <= ENDCOL(srow) ; col++)
        {
            tcnt = cnt ;
            PIRO_BAND_FLOPS(6*(srow-terow+1)) ;
            for (row = srow ; row >= terow ; row--)
            {
                ASSIGN_CS_SCALARS(c, s, givens, tcnt) ;
                tcnt++ ;
                PRINT(("Apply rot="ID"\n", tcnt)) ;
                PRINT(("row="ID"-"ID" col="ID"\n", row-1, row, col)) ;
                ASSIGN_TO_SCALAR(da, A, INDEX(row-1, col)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(row, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(row-1, col)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(row, col)) ;
            }
            if (col >= ENDCOL(erow-1))
            {
                terow++ ;
            }
        }

        tcnt = cnt ; /* Note the cnt before changing it */
        cnt = cnt + (srow-erow)+1 ;
        if (ENDCOL(srow) == ENDCOL(erow-1))
        {
            /* No fill possible from this wave of rotations */
            erow++ ;
            srow = MIN (srow+1, m-1) ;
            continue ;
        }

        /*-------------------- Phase 2 ----------------- */
        /* Create and remove fill immediately using a column rotation. Apply
         * new column rotation to entire column.*/
        col = ENDCOL(srow) ;
        for (row = srow ; row >= erow ; row--)
        {
            if (ENDCOL(row-1) == ENDCOL(row))
            {
                tcnt++ ;
                continue ;
            }

            ASSIGN_CS_SCALARS(c, s, givens, tcnt) ;
            tcnt++ ;
            PRINT(("Apply rot="ID"\n", tcnt)) ;
            PRINT(("row="ID"-"ID" col="ID"\n", row-1, row, col)) ;

            ASSIGN_TO_ZERO(da) ;
            ASSIGN_TO_SCALAR(db, A, INDEX(row, col)) ;
            APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
            ASSIGN_TO_MATRIX(db, A, INDEX(row, col)) ;
            ASSIGN_TO_SCALAR(ws, da, 0) ;

            ASSIGN_TO_SCALAR(da, A, INDEX(row-1, col-1)) ;
            ASSIGN_TO_SCALAR(db, ws, 0) ;

            GIVENS(da, db, d, c, s, Common) ;
            ASSERT (cgcnt < tcnt) ;
            SAVE_CS_SCALARS(givens, cgcnt, c, s) ;
            cgcnt++ ;

            ASSIGN_TO_MATRIX(d, A, INDEX(row-1, col-1)) ;

            PRINT(("Apply crot="ID"\n", cgcnt)) ;
            PRINT(("row="ID"-"ID" col="ID"-"ID"\n", row-1, derow, col-1, col)) ;
            /* Apply the rotation to the rest of the column */
            PIRO_BAND_FLOPS(6*(derow-row+1)+12) ;
            for (crow = row ; crow <= derow ; crow++)
            {
                ASSIGN_TO_SCALAR(da, A, INDEX(crow, col-1)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(crow, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(crow, col-1)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(crow, col)) ;
            }

            if (V != NULL)
            {
                s1 = MIN (col-c_srow, upd->vt_srow_min) ;
                APPLY_GIVENS_TO_COL_M(V, ldv, col-1, col, s1, c_lrow, c, s,
                                    i, da, db, dtemp, da1, db1, dtemp1, n) ;
            }

            col-- ;
        }

        erow++ ;
        srow = MIN (srow+1, m-1) ;

        if (V != NULL)
        {
            c_lrow = MIN (c_lrow+1, n-1) ;
            c_srow++ ;
        }
    } /* F blk */
    ASSERT (cnt == rgcnt) ;
    /* return the number of column rotations in this block */
    return cgcnt ;
}

/* ========================================================================== */
/* === apply_dblk =========================================================== */
/* ========================================================================== */

/* The D block is defined by rows dsrow..derow and columns bscol..becol
 * Apply row rotations, and remove any fill w/ column rotations and apply those.
 *
 * For every wave(row/column) of rotations that comes into the D block we will
 * need two passes on the D block.
 * Pass one will not generate any fill, but will apply all the column rotations
 *  w/o any fill.
 * Pass two will create and remove the fill using a column rotation and will
 * apply postpone applying the column rotation immediately, until all fill from
 * this wave is removed. Finally we can apply all row rotations as column
 * operations.
 *
 * Algorithm :
 *  1. For each wave of column rotations
 *      a. for each column in D-block (right-left)
 *          Apply current column rotation from current wave to column, generate
 *          fillin in current column.
 *          find row rotation to remove fill,
 *          apply row rotation only to remove fill (not to entire row)
 *      b. for each column in D-block (left-right)
 *          Apply the pending row rotations to current column.
 *      c. Adjust column indices for the next wave.
 *
 * Returns the number of new row rotations in the block.
 * */
Int PIRO_BAND(apply_dblk)
(
    Int m,              /* #rows in the matrix              */
    Int n,              /* #columns in the matrix           */
    Int nrc,            /* #columns in C                    */
    Int scol,           /* start column of the first wave   */
    Int ecol,           /* end column of the first wave     */
    Int cgcnt,          /* #rotations to be applied         */
    Int dsrow,          /* start row of the D block         */
    Int becol,          /* end column of the D block        */
    Int bl,             /* lower bandwidth of the matrix    */
    Int obu,            /* original lower b/w of A, usage hidden : A() macro */
    Int ldab,           /* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix                      */
    Entry *givens,      /* workspace for the rotations      */
    Entry *U,           /* o/p accumulated left rotations          */
    Int ldu,            /* leading dimension of u                  */
    Entry *C,           /* o/p for left rotations                  */
    Int ldc,            /* leading dimension of C                  */
    Int sym,            /* Is the matrix symmetric ?               */
    PIRO_BAND_LONG_NAME(uv_update) *upd, /* Structure for U and VT update */
    PIRO_BAND_Common *Common             /* common data structure         */
)
{
    Int rgcnt ;
    Int terow, terow1, terow2 ;
    Int row, col ;
    Int cnt, tcnt, rcnt ;
    Int i ;
    /* temproary variable for uv update */
    Int c_lrow = 0 ;
    Int c_srow = 0 ;
    Int s1 ;
    Entry c[CSIZE(1)], s[CSIZE(1)] ;
    Entry dtemp[CSIZE(1)], ws[CSIZE(1)] ;
    Entry d[CSIZE(1)] ;
    Entry da[CSIZE(1)], db[CSIZE(1)] ;
    Entry conj_s[2] ;
    Entry da1[CSIZE(1)], db1[CSIZE(1)], dtemp1[CSIZE(1)] ;
#ifdef PIROBAND_COMPLEX
    Entry conjs[CSIZE(1)], conjs1[CSIZE(1)] ;
#else
    Entry temp2[CSIZE(1)] ;
#endif

    PRINT(("D block\n")) ;
    rgcnt = 0 ;
    /* for U update */
    if (U != NULL)
    {
        c_lrow = upd->c_lrow ;
        c_srow = upd->k + 1 ;
    }
    ASSIGN_TO_ZERO(ws) ;

    for (cnt = 0 ; cnt < cgcnt ; )
    {
        /*-------------------- Phase 1 ----------------- */
        /* Apply the column rotation in current wave, generating any fill,
         * Remove the fill with a row rotation. */
        ASSERT (scol >= ecol && scol <= becol) ;
        PRINT(("scol="ID" ecol="ID" becol="ID"\n", scol, ecol, becol));
        for (col = scol ; col >= ecol ; col--)
        {
            terow1 = MIN (col+bl, m-1) ;
            terow2 = MIN (col-1+bl, m-1) ;
            ASSIGN_CS_SCALARS(c, s, givens, cnt) ;
            cnt++ ;
            PRINT(("Apply crot="ID"\n", cnt)) ;
            PRINT(("row="ID"-"ID" col="ID"-"ID"\n", dsrow, terow1, col-1, col));

            if (sym)
            {
                /* ASSIGN_TO_SCALAR(ws, A, INDEX(terow2, col)) ; */
                /* Using conj_s as a temp variable */
                ASSIGN_TO_SCALAR(conj_s, A, INDEX(terow2, col)) ;
                CONJ(ws, conj_s) ;
            }

            PIRO_BAND_FLOPS(6*(terow2-dsrow+1)) ;
            for (row = dsrow ; row <= terow2 ; row++)
            {
                ASSIGN_TO_SCALAR(da, A, INDEX(row, col-1)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(row, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(row, col-1)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(row, col)) ;
            }
            if (terow1 == terow2)
                continue ;
            ASSERT (terow2+1 == terow1) ;

            if (sym)
            {
                ASSIGN_TO_SCALAR(da, ws, 0) ;
            }
            else
            {
                ASSIGN_TO_ZERO(da) ;
            }
            ASSIGN_TO_SCALAR(db, A, INDEX(terow1, col)) ;
            APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
            ASSIGN_TO_SCALAR(ws, da, 0) ;
            ASSIGN_TO_MATRIX(db, A, INDEX(terow1, col)) ;

            ASSIGN_TO_SCALAR(da, A, INDEX(terow2, col-1)) ;
            ASSIGN_TO_SCALAR(db, ws, 0) ;

            PRINT(("Find rot="ID"\n", rgcnt)) ;
            PRINT(("row="ID"-"ID" col="ID"\n", terow2, terow1, col-1));

            if (sym)
            {
                /* Using conj_s as a temp variable */
                ASSIGN_TO_SCALAR(conj_s, s, 0) ;
                CONJ(s, conj_s) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(terow2, col-1)) ;

            }
            else
            {
                GIVENS(da, db, d, c, s, Common) ;
                ASSIGN_TO_MATRIX(d, A, INDEX(terow2, col-1)) ;
            }
            PIRO_BAND_FLOPS(12) ;

            SAVE_CS_SCALARS(givens, rgcnt, c, s) ;
            ASSERT (rgcnt < cnt) ;
            rgcnt++ ;

            CONJ(conj_s, s) ;
            if (U != NULL)
            {
                s1 = MIN (terow1-c_srow, upd->u_srow_min) ;
                APPLY_GIVENS_TO_COL_M(U, ldu, terow2, terow1, s1, c_lrow, c,
                    conj_s, i, da, db, dtemp, da1, db1, dtemp1, m) ;
            }
            if (C != NULL)
            {
                APPLY_GIVENS_TO_COL_M(C, ldc, terow2, terow1, 0, nrc-1, c,
                    conj_s, i, da, db, dtemp, da1, db1, dtemp1, m) ;
            }
        }

        if (U != NULL)
        {
            c_lrow = MIN (c_lrow+1, m-1) ;
            c_srow++ ;
        }

        if(MIN (scol+bl, m-1) == MIN (ecol-1+bl, m-1))
        {
            /* no fill possible in current wave */
            /* can use rgcnt > prevcnt, or fill=0/1 */
            scol = MIN (scol+1, n-1) ;
            ecol++ ;
            continue ;
        }

        /*-------------------- Phase 2 ----------------- */
        /* Apply the pending row rotations to the columns in this block.
         * */
        terow = MIN (ecol+bl, m-1) ;
        tcnt = rgcnt-1 ;
        for (col = ecol ; col <= becol ; col++)
        {
            row = terow ;
            ASSERT (tcnt >= 0 && tcnt < rgcnt) ;
            PIRO_BAND_FLOPS(6*(rgcnt-tcnt)) ;
            for (rcnt = tcnt ; rcnt < rgcnt ; rcnt++, row--)
            {
                ASSIGN_CS_SCALARS(c, s, givens, rcnt) ;
                PRINT(("Apply rot="ID"\n", rcnt)) ;
                PRINT(("row="ID"-"ID" col="ID"\n", row-1, row, col));

                ASSIGN_TO_SCALAR(da, A, INDEX(row-1, col)) ;
                ASSIGN_TO_SCALAR(db, A, INDEX(row, col)) ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ASSIGN_TO_MATRIX(da, A, INDEX(row-1, col)) ;
                ASSIGN_TO_MATRIX(db, A, INDEX(row, col)) ;
            }
            if (terow < m-1 && col < scol)
            {
                tcnt-- ;
                terow++ ;
            }
        }
        scol = MIN (scol+1, n-1) ;
        ecol++ ;

    }
    /* return the number of row rotations in this block */
    return rgcnt ;
}

/* ========================================================================== */
/* === set_to_eye =========================================================== */
/* ========================================================================== */

/* Initializa A, a mxn dense matrix, to identity matrix */
void PIRO_BAND(set_to_eye)
(
    Int m,          /* #rows in A */
    Int n,          /* #columns in A */
    Entry *A        /* A to set to identity */
)
{
    Int i, j ;
    for (j = 0 ; j < n ; j++)
    {
        for (i = 0 ; i < m ; i ++)
        {
            ASSIGN_ZERO_TO_MATRIX(A, CSIZE(i+j*m)) ;
        }
    }
    for (j = 0 ; j < MIN (m, n) ; j++)
    {
        A[CSIZE(j+j*m)] = ONE ;
    }
}

/* ========================================================================== */
/* === lowerband_transpose ================================================== */
/* ========================================================================== */

/* Simple band transpose of symmetric lower triangular band matrix to its
 * symmetric upper triangular band matrix. This is not a general band transpose. */
void PIRO_BAND(lowerband_transpose)
(
    Int n,          /* # columns in A           */
    Int bl,         /* lower bandwidth of A     */
    Entry *A,       /* A in band format         */
    Int lda,        /* leading dimension of A   */
    Entry *AT       /* A Transpose in band format */
)
{
    Int j ;
    Int Aindex, ATindex, cnt, entries ;
    Entry da[CSIZE(1)], ca[CSIZE(1)] ;

    ASSERT (bl >= 0 ) ;

    for (j = 0 ; j < n ; j++)
    {
        Aindex = j * lda ;
        ATindex = (j + 1) * (bl+1) - 1 ;
        entries = (j < n-bl) ? (bl+1) : (n-j) ;
        for (cnt = 0 ; cnt < entries ; cnt++)
        {
            /* Stride 1 access in A, Non stride 1 access in AT */
            ASSIGN_TO_SCALAR(da, A, CSIZE(Aindex)) ;
            CONJ(ca, da) ;
            ASSIGN_TO_MATRIX(ca, AT, CSIZE(ATindex)) ;
            Aindex++ ;
            ATindex += bl ;
        }
    }
}

/* ========================================================================== */
/* === general_transpose ==================================================== */
/* ========================================================================== */

/* Simple out of place dense transpose function for a dense matrix  */
void PIRO_BAND(general_transpose)
(
    Int m,          /* # rows in A              */
    Int n,          /* # columns in A           */
    Entry *A,       /* A to be transposed       */
    Int lda,        /* leading dimension of A   */
    Entry *AT,      /* A Transpose              */
    Int ldat        /* leading dimension of A'  */
)
{
    Int i, j ;
    Entry dtemp[2] ;
    Entry c_temp[2] ;

    for (j = 0 ; j < n ; j++)
    {
        for (i = 0 ; i < m ; i++)
        {
            ASSIGN_TO_SCALAR(dtemp, A, AINDEX(i, j, lda)) ;
            CONJ(c_temp, dtemp) ;
            ASSIGN_TO_MATRIX(c_temp, AT, AINDEX(j, i, ldat)) ;
        }
    }

}

/* ========================================================================== */
/* === inplace_conjugate_transpose ========================================== */
/* ========================================================================== */

/* Simple in place dense conjugate transpose for square matrices.  */
void PIRO_BAND(inplace_conjugate_transpose)
(
    Int n,          /* # columns in A           */
    Entry *A,       /* A to be transposed       */
    Int lda         /* leading dimension of A   */
)
{
    Int i, j ;
    Entry dtemp[2] ;
    Entry c_temp[2] ;

    for (j = 0 ; j < n ; j++)
    {
        ASSIGN_TO_SCALAR(dtemp, A, AINDEX(j, j, lda)) ;
        CONJ(c_temp, dtemp) ;
        ASSIGN_TO_MATRIX(c_temp, A, AINDEX(j, j, lda)) ;
        for (i = j+1 ; i < n ; i++)
        {
            ASSIGN_TO_SCALAR(dtemp, A, AINDEX(j, i, lda)) ;
            CONJ(dtemp, dtemp) ;
            ASSIGN_TO_SCALAR(c_temp, A, AINDEX(i, j, lda)) ;
            CONJ(c_temp, c_temp) ;
            ASSIGN_TO_MATRIX(dtemp, A, AINDEX(i, j, lda)) ;
            ASSIGN_TO_MATRIX(c_temp, A, AINDEX(j, i, lda)) ;
        }
    }
}

/* ========================================================================== */
/* === add_bidiag =========================================================== */
/* ========================================================================== */

/* Add the upper bidiagonal to a banded matrix A stored in packed band format.
 * LAPACK allows input matrices with no upper diagonal. This function adds an
 * upper bidiagonal to the input matrix A with lower bandwidth of kl. */
void PIRO_BAND(add_bidiag)
(
    Int n,              /* #columns in A                                 */
    Entry *A,           /* input matrix A                                */
    Int lda,            /* leading dimension of A                        */
    Int kl,             /* lower bandwidth of A                          */
    Entry *newA         /* output matrix with one upper bidiagonal added */
)
{
    Int i, j ;
    for (j = 0 ; j < n ; j++)
    {
        ASSIGN_ZERO_TO_MATRIX(newA, CSIZE(j*(kl+2))) ;
        for (i = 1 ; i < kl+2 ; i++)
        {
            ASSIGN_MATRIX_TO_MATRIX(newA, CSIZE(j*(kl+2)+i),
                                       A, CSIZE(j*lda+i-1)) ;
        }
    }
}

#ifdef PIROBAND_COMPLEX

/* ========================================================================== */
/* === complex_to_real ====================================================== */
/* ========================================================================== */

/* Converts the complex bidiagonal values to real and apply the corresponding
 * transformation to U, V, and C.
 * */

void PIRO_BAND(complex_to_real)
(
    Entry *A,               /* Band Matrix                             */
    Int ldab,               /* leading dimension of A                  */
    Int m,                  /* #rows in the original matrix            */
    Int n,                  /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,                 /* lower bandwidth                         */
    Int bu,                 /* upper bandwidth                         */
    Entry *U,               /* o/p accumulated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumulated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym                 /* Is the matrix symmetric ?               */
)
{
    Entry T[2] ;            /* temp variable       */
    Entry abst ;            /* absolute value of T */
    Entry temp[2], conj_T[2], res[2] ;
    Int obu ;
    Int i, j ;

    obu = bu ;

    if (!sym)
    {
        ASSIGN_TO_SCALAR(T, A, INDEX(0, 0)) ;
        for (i = 0 ; i < (MIN (m, n)) ; i++)
        {
            /* Assign the absolute value to the diagonal */
            abst = PIRO_BAND(hypot)(T[0], T[1]) ;
            A[INDEX(i, i)] = abst ;
            A[INDEX(i, i)+1] = 0.0 ;
            if (abst != 0.0)
            {
                T[0] = T[0]/abst ;
                T[1] = T[1]/abst ;
            }
            else
            {
                T[0] = 1 ;
                T[1] = 0.0 ;
            }

            /* U (1:m, i) = T * U(1:m, i) ; */
            if (U != NULL)
            {
                for (j = 0 ; j < m ; j++)
                {
                    ASSIGN_TO_SCALAR(temp, U, AINDEX(j, i, ldu)) ;
                    MULT(res, temp, T) ;
                    ASSIGN_TO_MATRIX(res, U, AINDEX(j, i, ldu)) ;
                }
            }

            if (C != NULL)
            {
                for (j = 0 ; j < nrc ; j++)
                {
                    ASSIGN_TO_SCALAR(temp, C, AINDEX(j, i, ldc)) ;
                    MULT(res, temp, T) ;
                    ASSIGN_TO_MATRIX(res, C, AINDEX(j, i, ldc)) ;
                }
            }

            if (i < (MIN (m, n))-1)
            {
                ASSERT (bu != 0) ;
                CONJ(conj_T, T) ;
                ASSIGN_TO_SCALAR(temp, A, INDEX(i, i+1)) ;
                MULT(T, temp, conj_T) ;
                abst = PIRO_BAND(hypot)(T[0], T[1]) ;
                /* Assign offdiagonal to absolute value of A(i, i+1) * T' */
                A[INDEX(i, i+1)] = abst ;
                A[INDEX(i, i+1)+1] = 0.0 ;
                if (abst != 0.0)
                {
                    T[0] = T[0]/abst ;
                    T[1] = T[1]/abst ;
                }
                else
                {
                    T[0] = 1 ;
                    T[1] = 0.0 ;
                }

                /* RT (1:n, i+1) = T' * RT (1:n, i+1) ; */
                CONJ(conj_T, T) ;
                if (V != NULL)
                {
                    for (j = 0; j < n ; j++)
                    {
                        ASSIGN_TO_SCALAR(temp, V, AINDEX(j, i+1, ldv)) ;
                        MULT(res, temp, conj_T) ;
                        ASSIGN_TO_MATRIX(res, V, AINDEX(j, i+1, ldv)) ;
                    }
                }
                /* T = A(i+1, i+1) * T' */
                ASSIGN_TO_SCALAR(temp, A, INDEX(i+1, i+1)) ;
                MULT(T, temp, conj_T) ;
            }
        }
    }
    else
    {
        /* Hermitian Matrix, Only the off diagonal elements are complex */
        for ( i = 1 ; i < n  ; i++)
        {
            /* Assign the absolute value to the off diagonal */
            ASSIGN_TO_SCALAR(T, A, INDEX(i-1, i)) ;
            abst = PIRO_BAND(hypot)(T[0], T[1]) ;
            A[INDEX(i-1, i)] = abst ;
            A[INDEX(i-1, i)+1] = 0.0 ;
            if (abst != 0.0)
            {
                T[0] = T[0]/abst ;
                T[1] = T[1]/abst ;
            }
            else
            {
                T[0] = 1 ;
                T[1] = 0.0 ;
            }
            if (i < n-1)
            {
                /* Adjust the next entry in the off diagonal */
                ASSIGN_TO_SCALAR(temp, A, INDEX(i, i+1)) ;
                MULT(res, temp, T) ;
                ASSIGN_TO_MATRIX(res, A, INDEX(i, i+1)) ;
            }

            CONJ(conj_T, T) ;
            if (U != NULL)
            {
                for (j = 0 ; j < m ; j++)
                {
                    ASSIGN_TO_SCALAR(temp, U, AINDEX(j, i, ldu)) ;
                    MULT(res, temp, conj_T) ;
                    ASSIGN_TO_MATRIX(res, U, AINDEX(j, i, ldu)) ;
                }
            }

            if (C != NULL)
            {
                for (j = 0 ; j < nrc ; j++)
                {
                    ASSIGN_TO_SCALAR(temp, C, AINDEX(j, i, ldc)) ;
                    MULT(res, temp, conj_T) ;
                    ASSIGN_TO_MATRIX(res, C, AINDEX(j, i, ldc)) ;
                }
            }

            if (V != NULL)
            {
                for (j = 0; j < n ; j++)
                {
                    ASSIGN_TO_SCALAR(temp, V, AINDEX(j, i, ldv)) ;
                    MULT(res, temp, conj_T) ;
                    ASSIGN_TO_MATRIX(res, V, AINDEX(j, i, ldv)) ;
                }
            }
        }
    }
}
#endif


/* ========================================================================== */
/* === piro_band_check ====================================================== */
/* ========================================================================== */

/* check if inputs are valid */

int PIRO_BAND(check)    /* return 0 if OK, < 0 if error */
(
    Int m,          /* number of rows in A */
    Int n,          /* number of columns in A */
    Int nrc,        /* number of rows in C */
    Int bl,         /* lower bandwidth */
    Int bu,         /* upper bandwidth */
    Entry *A,       /* the band matrix to be reduced */
    Int ldab,       /* leading dimension of A */
    Entry *B1,      /* the diagonal of the bidiagonal matrix */
    Entry *B2,      /* the superdiagonal of the bidiagional matrix */
    Entry *U,       /* accumulated left rotations */
    Int ldu,        /* leading dimension of U */
    Entry *V,       /* accumulated right rotations */
    Int ldv,        /* leading dimension of V */
    Entry *C,       /* for apply left rotations */
    Int ldc,        /* leading dimension of C */
    Int sym,        /* nonzero if A is symmetric, zero if A is unsymmetric */
    Int ccode       /* 0: LAPACK interface, nonzero: C interface */
)
{
    if (m == 0 || n == 0)
    {
        return (PIRO_BAND_OK) ;
    }
    if (m < 0)
    {
        return (PIRO_BAND_M_INVALID) ;
    }
    if (n < 0)
    {
        return (PIRO_BAND_N_INVALID) ;
    }
    if (sym && m != n)
    {
        return (PIRO_BAND_SYM_INVALID) ;
    }
    if (nrc < 0 || (C != NULL && nrc <= 0))
    {
        return (PIRO_BAND_NRC_INVALID) ;
    }
    if (bl < 0 || (sym && bl > 0) || bl > m-1)
    {
        return (PIRO_BAND_BL_INVALID) ;
    }
    if (bu < (ccode ? 1 : 0) || bu > n-1)
    {
        /* C interface requires bu > 0.  LAPACK requires bu >= 0 */
        return (PIRO_BAND_BU_INVALID) ;
    }
    if (!sym && ldab < bl+bu+1)
    {
        return (PIRO_BAND_LDAB_INVALID) ;
    }
    if (ldu < 0 || ( U != NULL && ldu < MAX (1, m)))
    {
        return (PIRO_BAND_LDU_INVALID) ;
    }
    if (ldv < 0 || ( V != NULL && ldv < MAX (1, n)))
    {
        return (PIRO_BAND_LDV_INVALID) ;
    }
    if (ldc < 0 || ( C != NULL && ldc < nrc))
    {
        return (PIRO_BAND_LDC_INVALID) ;
    }
    if (A == NULL)
    {
        return (PIRO_BAND_A_INVALID) ;
    }
    if (B1 == NULL)
    {
        return (PIRO_BAND_B1_INVALID) ;
    }
    if (B2 == NULL)
    {
        return (PIRO_BAND_B2_INVALID) ;
    }
    return (PIRO_BAND_OK) ;
}
