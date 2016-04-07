#include <stdlib.h>
#include "genband_util.h"
#include "genband_internal.h"
#include "blksky.h"

/* TODO:  i must always be a row index.  j must always be a column index.
   k can be either a row or column index.  Use p for a value from Ap or cp.
   Fix variable names accordingly. */

/* ========================================================================== */

/* Determine the top row of each column of A, and check the input matrix.
   On output, trow[j] is the top row index in column j.  Entries in the
   diagonal and +1 diagonal of A might not be in the data structure for A,
   but they are assumed to be present in the computation of trow[j].  trow[j]
   is also modified to ensure that entries A(i,j) where 0 <= j-i <= bw are
   present in the skyline data structure, even if they are not present in A.
   The row indices of A need not appear in sorted order.
   */

static Int init_col_list    /* returns 0 on success, -1 on error */
(
    /* input */
    sky_common *scom,
    Int *Ap,            /* size n+1, column pointers of A */
    Int *Ai,            /* size nz = Ap[n], row indices of A */
    Int n,

    /* output */
    Int *next,          /* size n+1 */
    Int *prev,          /* size n+1 */
    Int *trow,          /* size n */
    Int *encol,         /* size n */
    Int *size           /* scalar */
)
{
    Int cumsize, csize, clen, j, j2, rindex, bw, imin, imax, p, i ;

    trow [0] = 0 ;
    cumsize = 1 ;
    bw = MAX (1, scom->bw_est) ;

    for (j = 0 ; j < n ; j++)
    {
        /* pretend the bidiagonal entries are all present */
        imin = MAX (j-1, 0) ;
        imax = j ;

        /* Get the first and last entries in the jth column of A.  Assume row
           indices might not be sorted. */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            imin = MIN (i, imin) ;
            imax = MAX (i, imax) ;
        }

        if (imin < 0 || imax > j)
        {
            /* error return: matrix is not upper triangular, or it has
               row indices that are out of range.  Note that this test does
               not detect duplicate row indices.  TODO: decide what to do if
               the matrix has duplicate row indices. */
            return (-1) ;
        }

        if (j > 0)
        {
            /* assume A (rindex:j,j) is present in the skyline */
            rindex = MAX (j-bw, 0) ;
            trow [j] = MIN (rindex, imin) ;
            clen = j - trow [j] + 1 ;
            /* Heurstic : 5 more entries per column */
            /* TODO: "5" should be a parameter in scom */
            /* TODO: "size" is ignored by the caller.  why? */
            csize = clen + 5 ;
            cumsize += csize ;
            if (cumsize < 0)
            {
                /* error return: integer overflow */
                return (-1) ;
            }
            for (i = trow [j] ; i <= j ; i++)
            {
                encol [i] = j ; /* TODO Is there a better way to do this */
            }
        }

        next [j] = j+1 ;
        prev [j] = j-1 ;
    }

    next[n-1] = -1 ;
    *size = cumsize ;

    return (0) ;        /* success */
}

/* ========================================================================== */

/* TODO: return ncorners instead of void */

static void find_all_corners    /* TODO: why is this duplicated in blksky.c? */
(
    /* input */
    sky_common *scom,
    Int n,              /* TODO: describe all inputs/outputs, like this */
    Int *trow,

    /* output */
    Int *list,          /* size WHAT? */
    Int *ncorners       /* scalar; number of corners */
)
{
    Int k ;
    Int nlist ;
    Int bw ;

    bw = scom->bw_est ;

    /* TBD : Use iscorner */
    nlist = 0 ;
    for (k = 2 ; k < n-1 ; k++)
    {
        if (trow[k-1] >= trow[k] && trow[k] < trow[k+1] && k - trow[k] > bw)
        {
            /* column k is a corner */
            list [nlist++] = k ;
        }
    }
    if (trow[n-2] >= trow[n-1] && (n-1 - trow[n-1]) > bw)
    {
        /* column n-1 is a corner */
        list [nlist++] = n-1 ;
    }
    *ncorners = nlist ;
}

/* ========================================================================== */

static Int iscorner             /* TODO: why is this duplicated in blksky.c? */
(
    sky_common *scom,
    Int *trow,
    Int col,
    Int n
)
{
    Int bw ;

    bw = scom->bw_est ;

    if (col <= 1) return 0 ;

    if (col == n-1)
    {
        return (trow[col-1] >= trow[col] && (col - trow[col]) > bw ) ;
    }
    else
    {
        return (trow[col-1] >= trow[col]
             && trow[col] < trow[col+1]
             && (col - trow[col]) > bw ) ;
    }
}

#if 0
#ifndef NPRINT
static void print_int_array
(
    Int n,
    Int *Arr,
    char *str
)
{
    Int i ;
    PRINT (("%s ********\n", str )) ;
    for (i = 0 ;i < n ; i++)
    {
        PRINT ((" %ld ", Arr[i])) ;
        if (i != 0 && i % 10 == 0) PRINT (("\n")) ;
    }
    PRINT (("\n")) ;
}
#endif
#endif

/* ========================================================================== */

void blksky_find_symbolic       /* TODO return error code */
(
    sky_common *scom,
    Int *Ap,
    Int *Ai,
    double *Ax,
    Int m,
    Int n,
    sky_symbolic *sym
) 
{
    Int *intmem ;
    Int *tmp ;
    Int *next, *prev, *trow, *encol ;
    Int *list, *list2, *list3 ;   
    Int *r1, *c1 ;
    Int *rnew, *cnew ;
    Int *min_trow ;
    Int nlist, nlist2, nlist3 ; 
    Int head, tail ;
    Int size ;
    Int i, j, k ;
    Int top3 ;       /* TBD : used for list3 in k-2 corner check, required ?? */
    Int erow, ecol, e, tk, top ;
    Int srow, maxecol, rots, olderow, crow ;
    Int rrot_count, rot_count ;
    Int r1size, c1size ;
    Int rindex ;
    Int col, fcol, ccol, lcol, lcol1, row ;
    Int nexterow, nextsrow, nextecol ;
    Int temp ;
    Int kk ;
    Int bw ;
    /*
    double ws ;
    double da, db, d, c, s, dtemp ;
    */
    double swap_cnt, giv_cnt, fl ;
    Int status ;

    /*
    printf ("size of Int %g size of ptr %g\n",
        (double) sizeof (Int), (double) sizeof (Int *)) ;
    PRINT(("Entering blksky_reduce \n")) ;
    */

    intmem = (Int *) MALLOC (((5 * n ) +1) * sizeof(Int)) ;

    /* TODO: return if intmem is NULL */

    tmp = intmem ;

    /* Assigning the doubly linked list */
    prev = sym->prev ;
    next = sym->next ;
    min_trow = sym->min_trow ;

    /* Stacks of corners */
    list = tmp ;
    tmp += n ; 
    list2 = tmp ;
    tmp += n ; 
    list3 = tmp ;
    tmp += n ;
    trow = tmp ;
    tmp += (n + 1) ;
    encol = tmp ;
    tmp += n ;

    /*
    printf ("here I am 0b\n") ;
    */

    /* Initialize the integer data structures */
    status = init_col_list (scom, Ap, Ai, n, next, prev, trow, encol, &size) ;
    if (status < 0)
    {
        /* TODO return error code */
        printf ("Error!\n") ;
        FREE (intmem) ;
        return ;
    }

    /* TODO: the "size" output of init_col_list is ignored.  Why? */

    /*size += 3 * n ;*/
    /* cp[n] = size - 1 ;*/ /* not the size  but actual index */

    trow[n] = m ;   /* TODO why is this not done in init_col_list? */

    head = 0 ; /* TBD : Not red ? */
    tail = n - 1 ;

    /* Set counters to zero */
    swap_cnt = 0.0 ;
    giv_cnt = 0.0 ;
    fl = 0.0 ;
    sym->swap_cnt = 0.0 ;
    sym->giv_cnt = 0.0 ;
    sym->flp_cnt = 0.0 ;

    /*
    PRINT_INT_ARRAY (n+1, trow, "trow") ;
    PRINT_INT_ARRAY (n, encol, "encol") ;
    PRINT_INT_ARRAY (n, next, "next") ;
    PRINT_INT_ARRAY (n, prev, "prev") ;
    PRINT_INT_ARRAY (n+1, cp, "cp") ;
    PRINT_INT_ARRAY (nlist, list, "list") ;
    printf("size = %d \n", size) ;
    printf("nlist = %d \n", nlist) ;
    printf ("here I am 0c\n") ;
    */

    bw = scom->bw_est ;

    /*
    PRINT_INT_ARRAY (n+1, trow, "trow") ;
    PRINT_INT_ARRAY (n, encol, "encol") ;
    */

    /* Initialize the lists and top */
    for (i = 0 ; i < n ; i++)
    {
        list2[i] = -1 ;
        list3[i] = -1 ;

        if (min_trow != NULL) 
        {
            min_trow[i] = trow[i] ;
        }

        /* TODO: why are trow and encol being copied?  They should be computed
           in place, not in temp workspace and then just copied */
        sym->trow[i] = trow[i] ;

        if (trow[i] != 0)
        {
            PRINT(("i= "ID" trow[i]="ID" bw = "ID"\n", i, trow[i], bw)) ;
            ASSERT(i - trow[i] >= bw, "trow[i] wrong") ;
        }

        sym->encol[i] = encol[i] ;

        PRINT(("i = "ID" encol[i]="ID" bw = "ID"\n", i, encol[i], bw)) ;
        if (encol[i] != n-1)
        {
            ASSERT(encol[i]-i >= bw , "encol[i] wrong") ;
        }
    }

    sym->trow[n] =  trow[n] ;

    /* Quick return */
    if (min_trow == NULL)
    {
        /* TODO: what is this doing here? */
        FREE (intmem) ;
        return ;
    }

    /* Givens row and column numbers */
    r1 = (Int *) MALLOC (n * sizeof(Int)) ;
    c1 = (Int *) MALLOC (n * sizeof(Int)) ;
    r1size = n ;
    c1size = n ;

    find_all_corners(scom, n, trow, list, &nlist) ;

    nlist2 = 0 ;
    nlist3 = 0 ;
    top3 = 0 ; 

    while (nlist > 0)
    {
        /* peek the corner */
        /*printf("Finding a new block \n") ;*/
        k = list[nlist-1] ;
        i = trow[k] ;

        ASSERT(k >= 0 , "k < 0 \n") ;
        ASSERT(i >= 0 , "i < 0 \n") ;
        /* Find the block */
        erow = i ;
        ecol = k ;
        tk = k ;

        /*
        printf("nlist = %d \n", nlist) ;
        PRINT_INT_ARRAY (nlist, list, "list") ;
        PRINT_INT_ARRAY (n+1, trow, "trow") ;
        PRINT_INT_ARRAY (n, encol, "encol") ;
        */

        /* columns in blk limited by maxecol and e */
        if (i == 0)
        {
            e = 0 ;
        }
        else
        {
            e = MAX(encol[i-1], i) ;
        }

        /* sweep to the left */
        srow = k-1 ;
        maxecol = encol[srow] ;
        erow = MAX(trow[tk-1], erow) ;
        rots = 1 ; /* TBD : chk */
        while ( srow > e+1 && erow < srow - 1 && encol[tk-2] > ecol)
        {
            if (trow[tk-1] > erow)
            {
                if (trow[tk-1] < srow)
                {
                    olderow = erow ;
                    erow = MAX(trow[tk-1], erow) ;
                    if (erow == srow-1)
                    {
                        erow = olderow ;
                        break ;
                    }
                }
                else
                {
                    break ;
                }
            }
            tk = tk - 1 ;
            srow = tk - 1 ;
            maxecol = encol[srow] ;
            rots++ ;
        }

        /* sweep to the right */
        tk = k ;
        crow = i ;
        while (crow < erow && rots > 1)
        {
            if (trow[tk+1] == crow+1)
            {
                if (ecol+1 < maxecol)
                {
                    ecol++ ;
                    tk++ ;
                    crow++ ;
                    continue ;
                }
                else
                {
                    break ;
                }
            }
            rots-- ;
            crow++ ;
        }

        /*printf("Sweep to rows below \n") ;*/
        /* sweep to the rows below trow(k) */
        while (erow < srow-1 && rots > 1 && trow[tk+1] >= erow+1)
        {
            if (trow[tk+1] == erow+1)
            {
                if (ecol+1 < maxecol)
                {
                    ecol++ ;
                    tk++ ;
                    erow++ ;
                    continue ;
                }
                else
                {
                    break ;
                }
            }
            erow++ ;
            rots-- ;
        }

        /* Block size is fixed for the current iteration */

        /*
        printf("erow = %d\n", erow) ;
        printf("srow = %d\n", srow) ;
        printf("ecol = %d\n", ecol) ;
        printf("maxecol = %d\n", maxecol) ;
        printf("tk = %d\n", tk) ;
        printf("e = %d\n", e) ;
        */

        nlist2 = 0 ;
        top = list[nlist-1] ;
        ASSERT(top >= 0 , "top < 0 \n") ;
        while (top > srow && top <= ecol)
        {
            list[nlist-1] = -1 ;
            nlist-- ;

            if (trow[top] <= erow)
            {
                list2[nlist2++] = top ;
            }
            else
            {
                ASSERT(0, "trow[top] > erow \n") ;
                return ;
            }
            if (nlist > 0)
            {
                top = list[nlist-1] ;
                ASSERT(top >= 0 , "top < 0 \n") ;
            }
            else
            {
                break ;
            }
        }

        /*
        printf("nlist2 = %d \n", nlist2) ;
        PRINT_INT_ARRAY (nlist2, list2, "list2") ;
        printf("nlist = %d \n", nlist) ;
        PRINT_INT_ARRAY (nlist, list, "list") ;
        */

        /* zero all the corners in the current block */
        rrot_count = 0 ;
        rot_count = 0 ;
        while (nlist2 > 0)
        {
            k = list2[nlist2-1] ;
            list2[nlist2-1] = -1 ;
            nlist2-- ;
            i = trow[k] ;

            /*printf("current corner %d\n", k) ;
            printf("toprow of  corner %d\n", i) ;*/
            ASSERT ( k > 1 && i >= 0, "Illegal corner row/col") ;

            if (trow[k] == trow[k-1])
            {
                fl += k - trow[k] + 1 ;
                trow[k] = trow[k] + 1 ;
                /* Not reqd, trow[k] just increased.
                min_trow[k] = MIN(min_trow[k], trow[k]) ;*/
            }
            else
            {
                swap_cnt++ ;
                temp = trow[k-1] ;
                trow[k-1] = trow[k] ;
                trow[k] = temp ;
                /* Not reqd, swap should be inplace */
                min_trow[k-1] = MIN(min_trow[k-1], trow[k-1]) ;
                min_trow[k] = MIN(min_trow[k], trow[k]) ;
            }

            /* trow[k-1] amd trow[k] are swapped above. Use them correctly */
            /* TBD : Check */
            for (kk = trow[k-1] ; kk < MIN(trow[k], trow[k+1]) ; kk++)
            {
                if (encol[kk] == k)
                {
                    encol[kk] = k - 1 ;
                }
            }

            /* find the row rotation to zero ws */
            if (rrot_count == r1size)
            {
                rnew = (Int *) REALLOC (r1, r1size * 2 * sizeof(Int)) ;
                if (rnew != NULL)
                {
                    r1 = rnew ;
                    r1size = r1size *  2 ;
                }
            }
            r1[rrot_count++] = k ;
            fl += encol[k] - k + 2 ;

            /* Update the stack */
            /* If k-2 is no longer a corner remove it from the top of any of the
             * stacks */
            if (nlist2 > 0 && list2[nlist2-1] == k-2)
            {
                if (!iscorner(scom, trow, k-2, n))
                {
                    /*printf("Removing corner %d \n", k-2) ;*/
                    list2[nlist2-1] = -1 ;
                    nlist2-- ;
                }
            }
            else if (nlist3 > 0 && nlist3 > top3 && list3[top3] == k-2)
            {
                if (!iscorner(scom, trow, k-2, n))
                {
                    /*printf("Removing corner %d \n", k-2) ;*/
                    list3[top3] = -1 ;
                    top3++ ;
                    if (nlist3 == top3)
                    {
                        top3 = 0 ;
                        nlist3 = 0 ;
                    }
                }
            }
            else if (nlist > 0 && list[nlist-1] == k-2)
            {
                if (!iscorner(scom, trow, k-2, n))
                {
                    /*printf("Removing corner %d \n", k-2) ;*/
                    list[nlist-1] = -1 ;
                    nlist-- ;
                }
            }

            /* If k+1 is corner withing the block add it to list2 else if it is 
             * outside the block add it to list3. If list3 has smaller columns
             * in it pop them to list before adding k+1 to list3 
             * */
            if (iscorner(scom, trow, k+1, n))
            {
                if (trow[k+1] > erow || k+1 < srow || k+1 > ecol)
                {
                    /* Add the corner to list3 */
                    if (nlist3 > top3 && nlist3 > 0 && k+1 > list3[nlist3-1])
                    {
                        /* pop everything to list */
                        while (nlist3 > top3)
                        {
                    /*printf("Adding %d to right stack\n",list3[nlist3-1]);*/
                            list[nlist++] = list3[nlist3-1] ;
                            list3[--nlist3] = -1 ;
                        }
                        nlist3 = 0 ;
                        top3 = 0 ;
                    }
                    list3[nlist3++] = k+1 ;
                }
                else
                {
                    /*printf("Adding corner %d \n", k+1) ;*/
                    list2[nlist2++] = k+1 ;
                }
            }
            
            /* If k-1 is a corner within the block add it to list2. If it is 
             * outside the block add it to list3. If list2 has smaller columns
             * in it pop them to list before adding k-1 to list3 
             * */
            if (iscorner(scom, trow, k-1, n))
            {
                if (trow[k-1] > erow || k-1 <= srow || k-1 > ecol)
                {
                    /* Add the corner to list3 */
                    if (nlist3 > top3 && nlist3 > 0 && k-1 > list3[nlist3-1])
                    {
                        /* pop everything to list */
                        while (nlist3 > top3)
                        {
                    /*printf("Adding %d to right stack\n",list3[nlist3-1]);*/
                            list[nlist++] = list3[nlist3-1] ;
                            list3[--nlist3] = -1 ;
                        }
                        nlist3 = 0 ;
                        top3 = 0 ;
                    }
                    list3[nlist3++] = k-1 ;
                }
                else
                {
                    /*printf("Adding corner %d \n", k-1) ;*/
                    list2[nlist2++] = k-1 ;
                }
            }
        }
        /* Don't count the swaps, But count the row rotations for fill from
         * swaps */
        giv_cnt = giv_cnt + (2*rrot_count) ;

        /*printf("nlist = %d \n", nlist) ;
        PRINT_INT_ARRAY (nlist, list, "list") ;
        printf("nlist2 = %d \n", nlist2) ;
        PRINT_INT_ARRAY (nlist2, list2, "list2") ;
        printf("nlist3 = %d \n", nlist3) ;
        PRINT_INT_ARRAY (nlist3, list3, "list3") ;

        printf("top3=%d\n", top3) ;*/
        /* pop pending corners to rmc stack list */
        while (nlist3 > top3)
        {
            list[nlist++] = list3[nlist3-1] ;
            ASSERT (list3[nlist3-1] >= 0, "list3[nlist3-1] < 0") ;
            list3[--nlist3] = -1 ;
        }
        nlist3 = 0 ;
        top3 = 0 ;
        while (rrot_count > 0 || rot_count > 0)
        {
            /*printf("rrot_count = %d\n", rrot_count) ;
            printf("rot_count = %d\n", rot_count) ;
            PRINT_INT_ARRAY (rrot_count, r1, "r1") ;
            PRINT_INT_ARRAY (rot_count, c1, "c1") ;*/

            /* Complete column rotations in D block, Find row rotations and 
             * apply them in this block 
             * */
            for (rindex = 0 ; rindex < rot_count ; rindex++)
            {
                col = c1[rindex] ;
                if (rrot_count == r1size)
                {
                    rnew = (Int *) REALLOC (r1, r1size * 2 * sizeof(Int)) ;
                    if (rnew != NULL)
                    {
                        r1 = rnew ;
                        r1size *= 2 ;
                    }
                }
                r1[rrot_count++] = col ;
                fl += encol[col] - col + 2 ;
            }
            if (rot_count > 0)
            {
                giv_cnt = giv_cnt + rrot_count ;
            }

            fcol = ecol ;
            
            /* Row rotations only R block */
            for (ccol = ecol+1 ; ccol < encol[srow] ; ccol++)
            {
                fcol = ccol ;
            }
            fcol++ ;

            rot_count = 0 ;
            if (fcol > n-1)
            {
                rrot_count = 0 ;
                continue ;
            }

            nexterow = ecol ;
            nextecol = -1 ;
            nextsrow = -1 ;

            for (rindex = 0 ; rindex < rrot_count ; rindex++)
            {
                row = r1[rindex] ;
                lcol = encol[row] ;
                lcol1 = encol[row-1] ;
                /*printf("rindex = %d\n", rindex) ;
                printf("row = %d\n", row) ;
                printf("lcol = %d\n", lcol) ;
                printf("lcol1 = %d\n", lcol1) ;
                PRINT_INT_ARRAY (rrot_count, r1, "r1") ;*/
                if (lcol1 < lcol)
                {
                    encol[row-1] = lcol ;
                    trow[lcol] = row - 1 ;
                    if (iscorner(scom, trow, lcol, n))
                    {
                        if (nextecol == -1 && nextsrow == -1)
                        {
                            nextecol = lcol ;
                            nextsrow = lcol - 1 ;
                        }
                        nextecol = MAX(lcol, nextecol) ;
                        nextsrow = MIN(nextsrow, lcol-1) ;
                        if (rot_count == c1size)
                        {
                            cnew = (Int *) 
                                REALLOC (c1, c1size * 2 * sizeof(Int)) ;
                            if (cnew != NULL)
                            {
                                c1 = cnew ;
                                c1size *= 2 ;
                            }
                        }
                        c1[rot_count++] = lcol ;
                        fl += lcol - trow[lcol] + 1 ;

                        /* set trow and end encol to old values */
                        encol[row-1] = lcol1 ;
                        trow[lcol] = row ;
                        /* The fill is chasable, min_trow will remain the same*/
                    }
                    else
                    {
                        min_trow[lcol] = MIN(min_trow[lcol], trow[lcol]) ;
                    }
                }
            }
            giv_cnt = giv_cnt + rot_count ;
            rrot_count = 0 ;
            srow = nextsrow ;
            erow = nexterow ;
            ecol = nextecol ;
            /*printf("nextsrow = %d\n", srow) ;
            printf("nexterow = %d\n", erow) ;
            printf("nextecol = %d\n", ecol) ;
            printf("Completed the chase\n") ;*/
        }

    }

    /*
    printf("nlist = %d \n", nlist) ;
    PRINT_INT_ARRAY (nlist, list, "list") ;
    PRINT_INT_ARRAY (n+1, trow, "trow") ;
    PRINT_INT_ARRAY (n, encol, "encol") ;
    PRINT_INT_ARRAY (n, min_trow, "min_trow") ;
    */

    sym->swap_cnt = swap_cnt ;
    sym->giv_cnt = giv_cnt - swap_cnt ;
    sym->flp_cnt = fl ;

    FREE(intmem) ;
    FREE(r1) ;
    FREE(c1) ;

    /*PRINT(("Leaving blksky_reduce \n")) ;*/
}
