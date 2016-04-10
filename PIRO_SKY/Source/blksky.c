/*
 * Starting a TBD list here.
 * Mar 6th, 2009 : Need a mini pack when swapping columns and/or moving columns.
 * A column should never have more than m/2 space. In opt == 3 that may be the 
 * case when k-1 gets moved to the end, k will have k-1's space. When k gets 
 * moved we should ignore additional space and just move k.
 * */

#include <stdlib.h>
#include "genband_util.h"
#include "genband_internal.h"
#include "blksky.h"
#include "piro_band.h"

static void find_all_corners
(
    sky_common *scom,
    Int n,
    Int *trow,
    Int *list,
    Int *ncorners
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
            list[nlist++] = k ;
        }
    }
    if (trow[n-2] >= trow[n-1] && (n-1 - trow[n-1]) > bw)
    {
            list[nlist++] = n-1 ;
    }
    *ncorners = nlist ;
}

static Int iscorner
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

#ifndef NPRINT
static void print_int_array
(
    Int n,
    Int *Arr,
    char *str
)
{
    Int i ;
    PRINT(("%s ********\n", str )) ;
    for (i = 0 ;i < n ; i++)
    {
        PRINT((" %ld ", Arr[i])) ;
        if (i != 0 && i % 10 == 0) PRINT(("\n")) ;
    }
    PRINT(("\n")) ;
}

static void print_skyline
(
    sky_symbolic *sym,
    sky_skyline *sky 
)
{
    Int i, k ;
    Int m, n ;
    Int *cp, *trow ;
    double *A ;

    m = sky->m ;
    n = sky->n ;
    cp = sky->cp ;
    trow = sym->trow ;
    A = sky->A ;

    PRINT(("Skyline Matrix\n")) ;
    for (k = 0 ; k < n ; k++)
    {
        PRINT(("Column %ld", k)) ;
        for (i = trow[k] ; i <= k ; i++)
        {
            PRINT((" %0.4f", A[cp[k] - (k-i)])) ;
        }
        PRINT(("\n")) ;
        fflush(stdout) ;
    }
}
#endif

static void delete_column
(
    Int col,
    Int *next,
    Int *prev, 
    Int *head,
    Int *tail
)
{
    Int cnext, cprev ;

    cnext = next[col] ;
    cprev = prev[col] ;
    next[col] = -1 ;
    prev[col] = -1 ;

    if (cprev != -1)
    {
        next[cprev] = cnext ;
    }
    else
    {
        *head = cnext ;
    }

    if (cnext != -1)
    {
        prev[cnext] = cprev ;
    }
    else
    {
        *tail = cprev ;
    }
}

/* inserts col between prev[pos] and pos */
static void insert_column
(
    Int col,
    Int pos,
    Int *next,
    Int *prev, 
    Int *head,
    Int *tail
)
{
    Int cprev ;

    ASSERT(pos != -1, "pos == -1") ;
    cprev = prev[pos] ;

    if (cprev != -1)
    {
        next[cprev] = col ;
        prev[col] = cprev ;
    }
    else
    {
        *head = col ;
        prev[col] = -1 ;
    }

    prev[pos] = col ;
    next[col] = pos ;
}

static void add_column
(
    Int col,
    Int *next,
    Int *prev, 
    Int *head,
    Int *tail
)
{
    Int ctail ;
    ctail = *tail ;

    prev[col] = ctail ;
    next[col] = -1 ;

    if (ctail != -1)
    {
        next[ctail] = col ;
    }
    else
    {
        *head = col ;
    }

    *tail = col ;
}

#ifdef MATLAB_MEX_FILE
#ifndef NPRINT

static void print_skyline_svd
(
    sky_symbolic *sym,
    sky_skyline *sky,
    double *sgood
)
{
    Int i, k ;
    Int m, n ;
    Int *cp, *trow ;
    double *A ;
    mxArray *Smatlab, *pargout[2] ;
    double *Smatrix, *Svalues ;
    double maxerr, err ;

    m = sky->m ;
    n = sky->n ;
    cp = sky->cp ;
    trow = sym->trow ;
    A = sky->A ;

    Smatlab = mxCreateDoubleMatrix(n, n, mxREAL) ;
    Smatrix = mxGetPr (Smatlab) ;

    for (k = 0 ; k < n ; k++)
    {
        for (i = trow[k] ; i <= k ; i++)
        {
            Smatrix[i+k*n] = A[cp[k] - (k-i)] ;
        }
    }

    mexCallMATLAB (0, (mxArray **) NULL, 1, &Smatlab, "disp") ;
    /* get the SVD */
    mexCallMATLAB (1, pargout, 1, &Smatlab, "svd") ;

    Svalues = mxGetPr (pargout [0]) ;

    if (sgood != (double *) NULL)
    {
        maxerr = 0.0 ;
        for (k = 0 ; k < n ; k++)
        {
            err = ABS (Svalues [k] - sgood [k]) ;
            maxerr = MAX (maxerr, err) ;
        }
        PRINT (("maxerr %g  sgood[0] %g relerr %g\n", maxerr, sgood [0],
                    maxerr / sgood [0])) ;
        fflush(stdout) ;
    }

    mxDestroyArray(pargout[0]) ;
    mxDestroyArray(Smatlab) ;

}

#endif

static void find_sgood
(
    sky_symbolic *sym,
    sky_skyline *sky,
    double **sgood
)
{
    Int i, k ;
    Int m, n ;
    Int *cp, *trow ;
    double *A ;
    mxArray *Smatlab, *pargout[2] ;
    double *Smatrix, *Svalues ;

    PRINT (("\n\n------------------------------------------- find_sgood:\n")) ;

    m = sky->m ;
    n = sky->n ;
    cp = sky->cp ;
    trow = sym->trow ;
    A = sky->A ;

    Smatlab = mxCreateDoubleMatrix(n, n, mxREAL) ;
    Smatrix = mxGetPr (Smatlab) ;

    for (k = 0 ; k < n ; k++)
    {
        /*PRINT(("Column %ld", k)) ;*/
        for (i = trow[k] ; i <= k ; i++)
        {
            /*PRINT((" %0.4f", A[cp[k] - (k-i)])) ;*/
            Smatrix[i+k*n] = A[cp[k] - (k-i)] ;
        }
        /*PRINT(("\n")) ;
        fflush(stdout) ;*/
    }

    /* get the SVD */
    mexCallMATLAB (1, pargout, 1, &Smatlab, "svd") ;

    Svalues = mxGetPr (pargout [0]) ;

    PRINT(("Allocating sgood\n")) ;
    *sgood = (double *) MALLOC(n * sizeof(double)) ;
    for (i = 0 ; i < n ; i++)
    {
        (*sgood) [i] = Svalues [i] ;
    }

    mxDestroyArray(pargout[0]) ;
    mxDestroyArray(Smatlab) ;

}

#endif

/* ========================================================================== */
/* sky_reduce:  ... */
/* ========================================================================== */

/* describe me
    
    TODO change sky_* to skylinesvd_*
*/

Int sky_reduce      /* TODO why return Int? */
(
    sky_common *scom,
    sky_symbolic *sym,
    sky_skyline *sky,
    sky_factor *sfact
)
{
    Int m, n ;
    Int *intmem ;
    Int *tmp ;
    Int *next, *prev, *trow, *cp, *encol ;
    Int *list, *list2, *list3 ;   
    Int *r1, *c1 ;
    Int *rnew, *cnew ;
#if 0
    Int *min_trow ; /* not reqd */  /* TODO what is it? */
#endif
    Int nlist, nlist2, nlist3 ; 
    Int head, tail ;
    Int size ;
    Int i, k ;
    Int top3 ;       /* TBD : used for list3 in k-2 corner check, required ?? */
    Int erow, ecol, e, tk, top ;
    Int srow, maxecol, rots, olderow, crow ;
    Int rrot_count, rot_count ;
    Int r1size, c1size ;
    Int cpk, cpk1 ;
    Int rindex, cindex ;
    Int col, fcol, ccol, lcol, lcol1, row ;
    Int nexterow, nextsrow, nextecol ;
    Int rowindex ;
    Int tindex ;
    Int k1index ;
    Int ccindex ;
    Int temp ;
    Int kk ;
    Int save_srow, save_erow, save_ecol ;
    Int avail, required, avail2 ;
    Int opt, rcount ;
    Int ksize, k1size ;
    Int kroom, k1room ;
    Int elbow, basic_required ;
    Int ne, pr ;
    Int lcol_size ;
    Int realloc_cnt ;
    Int bw ;
    Int obu, ldab ;
    double *tempcol ;
    double ws ;
    double *A ;
    double *G1, *G2, *Gnew ;
    double da, db, d, c, s ;
    /*
    only used in variants of APPLY_GIVENS:
    double dtemp, conjs, temp2 ;
    double da1, db1, dtemp1 ;
    */
    Int tmpi ;
    double *sgood ;
    double *band_A ;

    /*
    printf ("size of Int %g size of ptr %g\n",
        (double) sizeof (Int), (double) sizeof (Int *)) ;
    printf ("sfact is %p\n", sfact) ;
    printf ("sfact->V is %p\n", sfact->V) ;
    if (sfact->V != NULL) printf ("sfact->V[0] is %g\n", sfact->V[0]) ;
    */

    /*PRINT(("Entering blksky_reduce \n"))) ;*/

    /* Initialize local pointers */
    opt = scom->opt ;
    rcount = scom->rcount ;
    scom->realloc_cnt = 0 ; 
    realloc_cnt = 0 ;

    trow = sym->trow ;
    encol = sym->encol ;
    prev = sym->prev ;
    next = sym->next ;
#if 0
    min_trow = sym->min_trow ;  /* TODO: never used */
#endif

    cp = sky->cp ;
    A = sky->A ;
    m = sky->m ;
    n = sky->n ;

    /* Allocate Integer memory for the stacks of corners */
    intmem = (Int *) MALLOC ((3 * n ) * sizeof(Int)) ;

    tmp = intmem ;

    /* Stacks of corners */
    list = tmp ;
    tmp += n ; 
    list2 = tmp ;
    tmp += n ; 
    list3 = tmp ;
    tmp += n ;/* not alloced */

    if (opt == 2)
    {
        /* Givens row and column numbers */
        r1 = (Int *) MALLOC (n * sizeof(Int)) ;
        c1 = (Int *) MALLOC (n * sizeof(Int)) ;
        r1size = n ;
        c1size = n ;

        /* Space to store the Givens rotations */
        G1 = (double *) MALLOC(2 * n * sizeof(double)) ;
        G2 = (double *) MALLOC(2 * n * sizeof(double)) ;
        rcount = n * n + 1 ;/* TBD : INT_MAX ?? , some large number actually */
    }
    else 
    {
        r1 = scom->iwork ;
        c1 = scom->iwork + rcount ;

        G1 = scom->dwork ;
        G2 = scom->dwork + (2 * rcount) ;
        r1size = rcount ;
        c1size = rcount ;
    }

    /* One temproary column, for the swap */
    tempcol = (double *) MALLOC(m * sizeof(double)) ;

    /* double space for the singular values */
    sgood = NULL ;

    size = sky->cp[n]+1 ;
    head = 0 ; /* TBD : Not reqd ?? */
    tail = n - 1 ;

    find_all_corners(scom, n, trow, list, &nlist) ;

#ifdef MATLAB_MEX_FILE
#ifndef NDEBUG
    find_sgood(sym, sky, &sgood) ;
#endif
#endif

    /*PRINT_SKYLINE (sym, sky) ;*/
    PRINT_INT_ARRAY (n+1, trow, "trow") ;
    PRINT_INT_ARRAY (n, encol, "encol") ;
    PRINT_INT_ARRAY (n, next, "next") ;
    PRINT_INT_ARRAY (n, prev, "prev") ;
    PRINT_INT_ARRAY (n+1, cp, "cp") ;
    PRINT_INT_ARRAY (nlist, list, "list") ;
    PRINT(("size = %ld \n", size)) ;
    PRINT(("nlist = %ld \n", nlist)) ;

    for (i = 0 ;i < n ; i++)
    {
        list2[i] = -1 ;
        list3[i] = -1 ;
        /*min_trow[i] = trow[i] ;*/
    }
    nlist2 = 0 ;
    nlist3 = 0 ;
    top3 = 0 ; 
    bw = scom->bw_est ;

    /* ---------------------------------------------------------------------- */
    /* reduce the matrix to bidiagonal form */
    /* ---------------------------------------------------------------------- */

    while (nlist > 0)       /* while there is a corner to eliminate ... */
    {

        PRINT (("\n------------------------------- skyline svd:\n")) ;
        PRINT_SKYLINE_SVD(sym, sky, sgood) ;
        PRINT (("\n\n\n")) ;

        /* peek the corner */
        PRINT(("Finding a new block \n")) ;
        k = list[nlist-1] ;
        i = trow[k] ;
        PRINT(("anchor column=%ld, its trow=%ld\n", k, i)) ;

        ASSERT(k >= 0 , "k < 0 \n") ;
        ASSERT(i >= 0 , "i < 0 \n") ;
        /* Find the block */
        erow = i ;
        ecol = k ;
        tk = k ;
        PRINT(("nlist = %ld \n", nlist)) ;
        PRINT_INT_ARRAY (nlist, list, "list") ;
        PRINT_INT_ARRAY (n+1, trow, "trow") ;
        PRINT_INT_ARRAY (n, encol, "encol") ;

        /* columns in blk limited by maxecol and e */
        if (i == 0)
        {
            e = 0 ;
        }
        else
        {
            e = MAX(encol[i-1], i) ;
        }


        /* ------------------------------------------------------------------ */
        /* sweep to the left */
        /* ------------------------------------------------------------------ */

        /*PRINT(("Sweep to left \n")) ;*/
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

        /* ------------------------------------------------------------------ */
        /* sweep to the right */
        /* ------------------------------------------------------------------ */

        /*PRINT(("Sweep to right \n")) ;*/
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

        /* ------------------------------------------------------------------ */
        /* sweep to the rows below trow(k) */
        /* ------------------------------------------------------------------ */

        /*PRINT(("Sweep to rows below \n")) ;*/
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

        /* ------------------------------------------------------------------ */
        /* Block size is now fixed for the current iteration */
        /* ------------------------------------------------------------------ */

        PRINT(("erow = %ld\n", erow)) ;
        PRINT(("srow = %ld\n", srow)) ;
        PRINT(("ecol = %ld\n", ecol)) ;
        PRINT(("maxecol = %ld\n", maxecol)) ;
        PRINT(("tk = %ld\n", tk)) ;
        PRINT(("e = %ld\n", e)) ;
        PRINT(("Block size is fixed \n")) ;
        PRINT_SKYLINE (sym, sky) ;

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
                return (0) ;
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
        PRINT(("nlist2 = %ld \n", nlist2)) ;
        PRINT_INT_ARRAY (nlist2, list2, "list2") ;
        PRINT(("nlist = %ld \n", nlist)) ;
        PRINT_INT_ARRAY (nlist, list, "list") ;

        /* ------------------------------------------------------------------ */
        /* zero all the corners in the current block */
        /* ------------------------------------------------------------------ */

        save_srow = srow ;
        save_erow = erow ;
        save_ecol = ecol ;
        while (nlist2 > 0)      /* TODO fix indenting */
        {
            rrot_count = 0 ;
            rot_count = 0 ;
            srow = save_srow ;
            erow = save_erow ;
            ecol = save_ecol ;

            while (nlist2 > 0 && rrot_count < rcount)
            {

            k = list2[nlist2-1] ;
            list2[nlist2-1] = -1 ;
            nlist2-- ;
            i = trow[k] ;
            PRINT(("Before corner \n")) ;
            PRINT_SKYLINE (sym, sky) ;

            PRINT(("current corner %ld\n", k)) ;
            PRINT(("toprow of  corner %ld\n", i)) ;
            ASSERT ( k > 1 && i >= 0, "Illegal corner row/col") ;

            if (trow[k] == trow[k-1])
            {
                da = A[cp[k-1] - ((k-1)-i)] ;
                db = A[cp[k] - (k-i)] ;
                ws = 0.0 ;
                PRINT(("Find & apply Givens r-%ld-%ld c-%ld-%ld",
                    i, k, k-1, k)) ;
                if (db != 0.0)
                {
                    GIVENS(da, db, d, c, s) ;
                    A[cp[k-1] - ((k-1)-i)]  = d ; /* HERE */
                    A[cp[k] - (k-i)] = 0.0 ;
                    cpk = cp[k] ;
                    cpk1 = cp[k-1] ;

                    for (rowindex = i+1 ; rowindex < k ; rowindex++)
                    {
                        da = A[cpk1 - ((k-1) - rowindex)] ;
                        db = A[cpk - (k - rowindex)] ;
                        APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                        A[cpk1 - ((k-1) - rowindex)] = da ;
                        A[cpk - (k - rowindex)]  = db ;
                    }
                    da = 0.0 ;
                    db = A[cpk] ;
                    APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                    ws = da ;
                    A[cpk] = db ;

                    if (sfact->V != NULL)
                    {
                        APPLY_GIVENS_TO_COL(sfact->V, n, k-1, k, 0, n-1, c, s, 
                            tmpi, da, db, dtemp, n) ;
                    }

                }
                trow[k] = trow[k] + 1 ;
                PRINT(("k= "ID" trow[k]="ID" bw = "ID"\n", k, trow[k], bw)) ;
                ASSERT(k - trow[k] >= bw, "trow[k] wrong") ;

            }
            else
            {

                PRINT(("swap rows %ld-%ld in columns %ld-%ld \n",
                    i, k, k, k-1)) ;
                if (opt != 3)
                {
                    cpk1 = cp[k-1] - ((k-1)-trow[k-1]) ;
                    for (rowindex = cpk1, tindex = 0 ; rowindex <= cp[k-1] ; 
                        rowindex++, tindex++)
                    {
                        tempcol[tindex] = A[rowindex] ;
                    }

                    cpk = cp[k] - (k-trow[k]) ;
                    k1index = cp[k-1] - ((k-1)-trow[k]) ;
                    for (ccindex = cpk ; ccindex < cp[k] ; ccindex++)
                    {
                        ASSERT(k1index > cp[k-2], "Trouble overwriting k-2") ;
                        A[k1index] = A[ccindex] ;
                        k1index++ ;
                    }

                    ws = A[cp[k]] ;

                    cpk = cp[k] - (k - trow[k-1]) ;
                    k1index = (k-1)-trow[k-1]+1 ;
                    for (ccindex = 0 ; ccindex < k1index ; ccindex++)
                    {
                        A[cpk] = tempcol[ccindex] ;
                        cpk++ ;
                    }
                    ASSERT(cpk == cp[k], "cpk != cp[k] while swap\n") ;
                    A[cp[k]] = 0.0 ;
                }
                else
                {
                    ksize = k - trow[k] + 1 ;
                    k1size = (k-1) - trow[k-1] + 1 ;
                    PRINT(("ksize = %ld, k1size = %ld \n", ksize, k1size)) ;

                    if (prev[k] == k-1)
                    {
                        /* k and and k-1 are adjacent in memory */

                        /* room on top of k (kroom) and k-1(k1room) */
                        kroom = cp[k] - cp[prev[k]] - ksize ;
                        k1room = cp[k-1] - cp[prev[k-1]] - k1size ;
                        PRINT(("kroom = %ld, k1room = %ld \n", kroom, k1room)) ;

                        /* NO OPTIMIZATIONS HERE , PLAIN SWAP */

                        /* cpk  - where we copy k to (index of new k-1).
                         * cpk1 - where we copy k-1 to (index of new k). */
                        cpk = cp[prev[k-1]]+1 ;
                        /* We copy ksize-1 entries to k-1, and diag to ws */
                        cpk1 = cpk + k1room + (ksize-1) ;
                        PRINT(("cpk = %ld, cpk1 = %ld \n", cpk, cpk1)) ;
                    }
                    else
                    {
                        /* Can k and k-1 be next to each other in memory ? */
                        basic_required = k1size + ksize ;

                        /* HERE assuming 10 more spaces again */
                        /* TBD : Should be based on max column height */
                        required = k1size + ksize + 10 ;
                        kroom = 5 ;
                        k1room = 5 ;

                        /* avail room near k-1 */
                        ne = next[k-1] ;
                        pr = prev[k-1] ;
                        if (ne == -1) 
                        {
                            /* k-1 is the tail. */
                            /* If k-1 is the tail we will not use the elbow
                             * room. (We can if we want in aggressive cases)
                             * We will try to squeeze k and k-1 in 
                             * and cp[prev[k-1]] .. cp[k-1] if possible.
                             * ************WON"t WORK NOW change delete_column 
                             * below
                             * */
                            ne = k-1 ;
                        }
                        ASSERT(ne >= 0 && pr >= 0 , "invalid ne/pr") ;
                        avail = cp[ne] - (ne-trow[ne]+1) - cp[pr] ;

                        ne = next[k] ;
                        pr = prev[k] ;
                        if (ne == -1) 
                        {
                            /* k is the tail. */
                            /* If k is the tail we will not use the elbow
                             * room. We will try to squeeze k-1 in between k
                             * and prev[k] if possible
                             * */
                            ne = k ;
                        }
                        ASSERT(ne >= 0 && pr >= 0 , "invalid ne/pr") ;
                        avail2 = cp[ne] - (ne-trow[ne]+1) - cp[pr] ;

                        /* How to split the available space ?
                         * 1. Equally between k and k-1.
                         * 2. 5 more space for k and k-1 each, rest to k-2
                         *  like a mini pack. (current) Will replace 5 
                         *  using max column ht. (TBD)
                         *  */
                        if (avail < required && avail2 < required)
                        {
                            /* Do each of them have space to atleast hold the
                             * other column ? If not
                             * need to use elbow room or realloc
                             * */
                            if (avail < basic_required && 
                                    avail2 < basic_required)
                            {
                                if (avail < ksize-1)
                                {
                                    /* space near k-1 cannot hold even entries
                                     * from k, move
                                     * k-1 to elbow or realloc. k will surely
                                     * hold k-1. */
                                    cpk1 = cp[prev[k]] + 1 ;
                                    kroom = (avail2 - k1size)/2 ;

                                    elbow = cp[n] - cp[tail] ;

                                    if (elbow >= ((ksize-1) + 5) )
                                    {
                                        k1room = 5 ;
                                        cpk = cp[tail] + 1 + k1room ;
                                    }
                                    else
                                    {
                                        /* If k-1 == tail we can use it but
                                         * later */
                                        Gnew = (double *) REALLOC(A, 
                                            ((cp[n]+1)+10*n)* sizeof(double)) ;
                                        if (Gnew != NULL)
                                        {
                                            A = Gnew ;
                                            sky->A = Gnew ;
                                        }
                                        else
                                        {
                                            PRINT(("Out of memory \n")) ;
                                            /* TBD :Should call pack here */
                                        }
                                        realloc_cnt++ ;
                                        for (tindex=size; tindex<size+10 *n ;
                                                tindex++)
                                        {
                                            A[tindex] = 0.0 ;
                                        }
                                        size = size + 10 * n ;
                                        /* Actual index not the size, 
                                         * cp[tail] != cp[n] */
                                        cp[n] = cp[n] + 10 * n ; /* HERE */
                                        k1room = 5 ;
                                        cpk = cp[tail] + 1 + k1room ;
                                    }

                                    /* Add k-1 to end of linked list */
                                    delete_column(k-1, next, prev, &head, 
                                            &tail) ;
                                    add_column(k-1, next, prev, &head, 
                                            &tail) ;
                                }
                                else
                                {
                            /* Leave them in different locations in memory */

                            /* This takes away half of whatever top room 
                             * next[k/k-1] had.  FIX */
                                    kroom = (avail2 - k1size)/2 ;
                                    k1room = (avail - (ksize-1))/2 ;
                                    cpk = cp[prev[k-1]]+ 1 ;
                                    cpk1 = cp[prev[k]] + 1 ;
                                }
                            }
                            else
                            {
                                /* Bring them together. Does this make sense. ?
                                 * Even if there is one more fill they will 
                                 * move away again if MAX(avail,avail2) ==
                                 * basic_required. May be the room should be 5
                                 * here when we we include max column height
                                 * above. room == 0 should be used only for
                                 * aggressive.
                                 * */
                                /* Choose the one with the maximum space */
                                pr = (avail > avail2 ? k-1 : k ) ;

                                /* Should use this after modifying 10/5 above
                                 * room can be 5 here rather than max col ht. */
                                /* tosplit = (avail > avail2 ? avail : avail2) ;
                                tosplit = tosplit - ksize - k1size ;*/

                                cpk = cp[prev[pr]]+1 ;
                                cpk1 = cpk + (ksize-1) ;
                                kroom = 0 ;
                                k1room = 0 ;

                                /* Need an insert column here in the order 
                                 * prev[pr] .. k-1 .. k .. next[pr].
                                 * */
                                if (pr == k-1)
                                {
                                    delete_column(k, next, prev, &head, &tail) ;
                                    if (k-1 == tail)
                                    {
                                        add_column(k, next, prev, &head, &tail);
                                    }
                                    else
                                    {
                                        insert_column(k, next[k-1], next, prev, 
                                                    &head, &tail);
                                    }
                                }
                                else
                                {
                                    delete_column(k-1, next, prev, &head, 
                                                &tail) ;
                                    insert_column(k-1, k, next, prev, &head, 
                                                &tail);
                                }

                            }
                        }
                        else 
                        {
                         /* Both locations have enough space to fit k and 
                          * k-1. Find the natural one */
                            /*if (prev[k-1] == k-2)
                            {
                                pr = k-1 ;
                            }
                            else if (next[k] == k+1)
                            {
                                pr = k ;
                            }
                            else
                            {*/
                                /* Choose the one with the maximum space */
                                pr = (avail > avail2 ? k-1 : k ) ;
                            /*}*/

                            /* cpk  - where we can copy k.
                             * cpk1 - where we can copy k-1. */
                            cpk = cp[prev[pr]]+1 ;
                            cpk1 = cpk + 5 + (ksize-1) ;

                            /* Need an insert column here in the order 
                             * prev[pr] .. k-1 .. k .. next[pr].
                             * */
                            if (pr == k-1)
                            {
                                delete_column(k, next, prev, &head, &tail) ;
                                if (k-1 == tail)
                                {
                                    add_column(k, next, prev, &head, &tail);
                                }
                                else
                                {
                                    insert_column(k, next[k-1], next, prev, 
                                                &head, &tail);
                                }
                            }
                            else
                            {
                                delete_column(k-1, next, prev, &head, &tail) ;
                                insert_column(k-1, k, next, prev, &head, &tail);
                            }

                        }
                    }

                    PRINT(("Copying %ld to tempcol\n", k-1)) ;
                    rowindex = cp[k-1] - k1size + 1 ;
                    for (tindex = 0 ; rowindex <= cp[k-1] ; )
                    {
                        ASSERT (tindex < m && rowindex < size, 
                                        "tindex/rowindex maxed") ;
                        ASSERT (tindex >= 0 && rowindex >= 0, 
                                        "tindex/rowindex < 0") ;
                        tempcol[tindex++] = A[rowindex++] ;
                    }

                    PRINT(("Copying %ld to %ld\n", k, k-1)) ;
                    ccindex = cp[k] - ksize + 1 ;
                    /* Need to add the right size to k1index */
                    for ( k1index = cpk + k1room ; ccindex < cp[k] ; )
                    {
                        PRINT(("k1index = %ld\n", k1index)) ;
                        ASSERT (k1index < size && ccindex < size,
                                        "k1index/ccindex maxed") ;
                        ASSERT (k1index >= 0 && ccindex >= 0,
                                        "k1index/ccindex < 0") ;
                        A[k1index++] = A[ccindex++] ;
                    }

                    cp[k-1] = k1index - 1 ;
                    ws = A[cp[k]] ;

                    /* Need to add the right size to k1index */
                    k1index = cpk1 + kroom ;
                    PRINT(("Copying tempcol to %ld\n",  k)) ;
                    for (ccindex = 0 ; ccindex < tindex ; ccindex++)
                    {
                        PRINT(("k1index = %ld\n", k1index)) ;
                        ASSERT (k1index < size && ccindex < m,
                                        "k1index/ccindex maxed") ;
                        ASSERT (k1index >= 0 && ccindex >= 0,
                                        "k1index/ccindex < 0") ;
                        A[k1index++] = tempcol[ccindex] ;
                    }
                    PRINT(("Copying 0 to diag\n")) ;
                    ASSERT (k1index < size && k1index >= 0, "k1index wrong") ;
                    A[k1index] = 0.0 ;

                    PRINT(("Setting right the cp array\n")) ;
                    cp[k] = k1index ;

                }

                PRINT(("Setting right the trow array\n")) ;

                temp = trow[k-1] ;
                trow[k-1] = trow[k] ;
                trow[k] = temp ;
                PRINT(("k= "ID" trow[k]="ID" bw = "ID"\n", k, trow[k], bw)) ;
                ASSERT(k - trow[k] >= bw, "trow[k] wrong") ;
                ASSERT((k-1) - trow[k-1] >= bw, "trow[k-1] wrong") ;

#if 0
                /* TODO min_trow is never used, except for assertions */
                if (min_trow != NULL)
                {
                ASSERT(trow[k-1] >= min_trow[k-1], "trow[k-1] unexpected") ;
                ASSERT(trow[k] >= min_trow[k], "trow[k] unexpected") ;
                }
#endif

                PRINT(("Did trow, now apply givens to col\n")) ;

                if (sfact->V != NULL)
                {
                    /*mexPrintf("I am here\n") ;*/
                    c = 0.0 ;
                    s = 1.0 ;

                    /*
                    printf ("sfact is %p\n", (void *) sfact) ;
                    printf ("sfact->V is %p\n", (void *) (sfact->V)) ;
                    PRINT(("Go to givens to col\n")) ;
                    */

                    APPLY_GIVENS_TO_COL(sfact->V, n, k-1, k, 
                        0, n-1, c, s, tmpi, da, db, dtemp, n) ;
                }

                PRINT(("Did apply givens to col\n")) ;

            }

            PRINT(("After corner \n")) ;
            PRINT_SKYLINE (sym, sky) ;

            PRINT(("head = %ld, tail=%ld\n", head, tail)) ;
            PRINT_INT_ARRAY (n, next, "next") ;
            PRINT_INT_ARRAY (n, prev, "prev") ;
            PRINT_INT_ARRAY (n+1, cp, "cp") ;
            PRINT(("tindex \n")) ;
            for (tindex = 0 ; tindex < n ; tindex++)
            {
                PRINT(("%ld ", tindex)) ;
#ifndef NDEBUG
                if (prev[tindex] != -1)
                    ASSERT(cp[tindex] > cp[prev[tindex]], "cp[ti]<=cp[p[ti]") ;
#endif
            }

            PRINT(("Reassigning encol values \n"));
            /* trow[k-1] amd trow[k] are swapped above. Use them correctly */
            /* TBD : Check */
            for (kk = trow[k-1] ; kk < MIN(trow[k], trow[k+1]) ; kk++)
            {
                if (encol[kk] == k)
                {
                    encol[kk] = k - 1 ;
                    PRINT(("kk= "ID" encol[kk]="ID" bw = "ID"\n",
                        kk, encol[kk], bw)) ;
                    if (encol[kk] != n-1)
                    {
                        ASSERT(encol[kk]-kk >= bw , "encol[k] wrong") ;
                    }
                }
            }

            PRINT(("Removing fill using row rotations \n"));

            if (ws != 0.0)
            {
                PRINT(("ws is %g\n", ws));

                da = A[cp[k-1]] ;
                db = ws ;
                d = 0 ;
                c = 0 ;
                s = 0 ;
                PRINT(("stuck? da %g db %g d %g c %g s %g\n", da, db, d, c, s));
                /* TBD : we can swap if da == 0 */
                GIVENS(da, db, d, c, s) ;
                PRINT(("did:   da %g db %g d %g c %g s %g\n", da, db, d, c, s));
                PRINT(("here 33\n"));
                A[cp[k-1]] = d ; /* HERE */
                ws = 0.0 ;

                PRINT(("Apply rowrow r-%ld-%ld c-%ld-%ld\n",
                    k-1, k, k-1, ecol)) ;
                PRINT(("rot no %ld\n", rrot_count)) ;

                for (cindex = k ; cindex <= ecol ; cindex++) 
                {
                    cpk = cp[cindex] ;
                    da = A[cpk - (cindex-(k-1))] ;
                    db = A[cpk - (cindex-k)] ;
                    APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                    A[cpk - (cindex-(k-1))] = da ;
                    A[cpk - (cindex-k)] = db ;
                }
                if (sfact->U != NULL)
                {
                    APPLY_GIVENS_TO_COL(sfact->U, m, k-1, k, 0, m-1, c, s, tmpi,
                        da, db, dtemp, m) ;
                }

                if (rrot_count == r1size)
                {
                    PRINT(("Realloc of r1 and G2 ***************")) ;
                    ASSERT(opt == 2, "opt != 2") ;
                    rnew = (Int *) REALLOC (r1, r1size * 2 * sizeof(Int)) ;
                    Gnew = (double *) REALLOC(G2, r1size*4*sizeof(double)) ;
                    if (rnew != NULL)
                    {
                        r1 = rnew ;
                        G2 = Gnew ;
                        r1size = r1size *  2 ;
                    }
                }
                G2[2*rrot_count] = c ;
                G2[2*rrot_count+1] = s ;
                r1[rrot_count++] = k ;
            }
            else
            {
                if (rrot_count == r1size)
                {
                    PRINT(("Realloc of r1 and G2 ***************")) ;
                    ASSERT(opt == 2, "opt != 2") ;
                    rnew = (Int *) REALLOC (r1, r1size * 2 * sizeof(Int)) ;
                    Gnew = (double *) REALLOC(G2, r1size*4*sizeof(double)) ;
                    if (rnew != NULL)
                    {
                        r1 = rnew ;
                        G2 = Gnew ;
                        r1size = r1size *  2 ;
                    }
                }
                G2[2*rrot_count] = 1 ;
                G2[2*rrot_count+1] = 0 ;
                r1[rrot_count++] = k ;
            }
            /* TBD : I am skipping the multiply by identity part in MATLAB */

            PRINT (("Update the stack\n")) ;

            /* Update the stack */
            /* If k-2 is no longer a corner remove it from the top of any of the
             * stacks */
            if (nlist2 > 0 && list2[nlist2-1] == k-2)
            {
                if (!iscorner(scom, trow, k-2, n))
                {
                    PRINT(("Removing corner %ld \n", k-2)) ;
                    list2[nlist2-1] = -1 ;
                    nlist2-- ;
                }
            }
            else if (nlist3 > 0 && nlist3 > top3 && list3[top3] == k-2)
            {
                if (!iscorner(scom, trow, k-2, n))
                {
                    PRINT(("Removing corner %ld \n", k-2)) ;
                    list3[top3] = -1 ;
                    top3++ ;
                    if (nlist3 == top3)
                    {
                        top3 = 0 ; /* HERE */
                        nlist3 = 0 ;
                    }
                }
            }
            else if (nlist > 0 && list[nlist-1] == k-2)
            {
                if (!iscorner(scom, trow, k-2, n))
                {
                    PRINT(("Removing corner %ld \n", k-2)) ;
                    list[nlist-1] = -1 ;
                    nlist-- ;
                }
            }

            PRINT (("yet more\n")) ;

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
                        while (nlist3 > top3) /* HERE : All top3 usage */
                        {
                            PRINT(("Adding %ld to right stack\n",
                                list3[nlist3-1]));
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
                    PRINT(("Adding corner %ld \n", k+1)) ;
                    list2[nlist2++] = k+1 ;
                }
            }

            PRINT (("yet more 2\n")) ;
            
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
                            PRINT(("Adding %ld to right stack\n",
                                list3[nlist3-1]));
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
                    PRINT(("Adding corner %ld \n", k-1)) ;
                    list2[nlist2++] = k-1 ;
                }
            }
        }

        /*PRINT(("nlist = %ld \n", nlist)) ;
        PRINT_INT_ARRAY (nlist, list, "list") ;
        PRINT(("nlist2 = %ld \n", nlist2)) ;
        PRINT_INT_ARRAY (nlist2, list2, "list2") ;
        PRINT(("nlist3 = %ld \n", nlist3)) ;
        PRINT_INT_ARRAY (nlist3, list3, "list3") ;*/

#if 0
        PRINT(("top3=%ld\n", top3)) ;
        /* pop pending corners to rmc stack list */
        while (nlist3 > top3)
        {
            PRINT(("Adding %ld to right stack\n", list3[nlist3-1])) ;
            list[nlist++] = list3[nlist3-1] ;
            ASSERT (list3[nlist3-1] >= 0, "list3[nlist3-1] < 0") ;
            list3[--nlist3] = -1 ;
        }
        nlist3 = 0 ;
        top3 = 0 ;
#endif
        PRINT(("Before the chase \n")) ;
        /*mexPrintf("rrot_count=%ld \n", rrot_count ) ;*/
        PRINT_SKYLINE (sym, sky) ;
        while (rrot_count > 0 || rot_count > 0)
        {
            PRINT(("rrot_count = %ld\n", rrot_count)) ;
            PRINT(("rot_count = %ld\n", rot_count)) ;
            PRINT_INT_ARRAY (rrot_count, r1, "r1") ;
            PRINT_INT_ARRAY (rot_count, c1, "c1") ;

            /* Column rotations only C block */
            if (erow+1 <= srow-1)
            {
                PRINT(("In Par col rotation block.\n")) ;
                for (rindex = 0 ; rindex < rot_count ; rindex++)
                {
                    c = G1[2*rindex] ;
                    s = G1[2*rindex+1] ;
                    col = c1[rindex] ;
                    PRINT(("Apply colrot r-%ld-%ld c-%ld-%ld\n",
                        erow+1, srow-1, col-1, col)) ;
                    PRINT(("rot no %ld\n", rindex)) ;

                    cpk = cp[col] ;
                    cpk1 = cp[col-1] ;
                    for (rowindex = erow+1 ; rowindex < srow ; rowindex++)
                    {
                        da = A[cpk1 - ((col-1)-rowindex)] ;
                        db = A[cpk - (col-rowindex)] ;
                        APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                        A[cpk1 - ((col-1)-rowindex)] = da ;
                        A[cpk - (col-rowindex)] = db ;
                    }
                }
            }



            PRINT(("In diag block\n")) ;
            /* Complete column rotations in D block, Find row rotations and 
             * apply them in this block 
             * */
            for (rindex = 0 ; rindex < rot_count ; rindex++)
            {
                c = G1[2*rindex] ;
                s = G1[2*rindex+1] ;
                col = c1[rindex] ;
                PRINT(("Apply colrot r-%ld-%ld c-%ld-%ld\n",
                    srow, col, col-1, col)) ;
                PRINT(("rot no %ld\n", rindex)) ;

                cpk = cp[col] ;
                cpk1 = cp[col-1] ;

                for (rowindex = srow ; rowindex < col ; rowindex++)
                {
                    da = A[cpk1 - ((col-1)-rowindex)] ;
                    db = A[cpk - (col-rowindex)] ;
                    APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                    A[cpk1 - ((col-1)-rowindex)] = da ;
                    A[cpk - (col-rowindex)] = db ;
                }

                /* Complete the column rotation by creating the fill */
                da = 0.0 ;
                db = A[cpk] ;
                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                ws = da ; 
                A[cpk] = db ;

                c1[rindex] = -1 ;/* TBD : may not be needed , Just in case */

                /* Find the row rotation to eliminate the fill */
                if (ws != 0.0)
                {
                    da = A[cpk1] ;
                    db = ws ;
                    GIVENS(da, db, d, c, s ) ;
                    A[cpk1] = d ; /* HERE */
                    ws = 0.0 ;

                    PRINT(("Apply rowrot r-%ld-%ld c-%ld-%ld\n",
                        col-1, col, col-1, ecol)) ;
                    PRINT(("rot no %ld\n", rrot_count)) ;

                    for (cindex = col ; cindex <= ecol ; cindex++) 
                    {
                        cpk = cp[cindex] ;
                        da = A[cpk - (cindex-(col-1))] ;
                        db = A[cpk - (cindex-col)] ;
                        APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                        A[cpk - (cindex-(col-1))] = da ;
                        A[cpk - (cindex-col)] = db ;
                    }

                    if (sfact->U != NULL)
                    {
                        APPLY_GIVENS_TO_COL(sfact->U, m, col-1, col, 0, m-1, c,
                            s, tmpi, da, db, dtemp, m) ;
                    }

                    if (rrot_count == r1size)
                    {
                        PRINT(("Realloc of r1 and G2 ***************")) ;
                        ASSERT(opt == 2, "opt != 2") ;
                        rnew = (Int *) REALLOC(r1, r1size * 2 * sizeof(Int)) ;
                        Gnew = (double *) REALLOC(G2, r1size*4*sizeof(double)) ;
                        if (rnew != NULL)
                        {
                            r1 = rnew ;
                            G2 = Gnew ;
                            r1size *= 2 ;
                        }
                    }
                    G2[2*rrot_count] = c ;
                    G2[2*rrot_count+1] = s ;
                    r1[rrot_count++] = col ;
                }
                else
                {
                    if (rrot_count == r1size)
                    {
                        PRINT(("Realloc of r1 and G2 ***************")) ;
                        ASSERT(opt == 2, "opt != 2") ;
                        rnew = (Int *) REALLOC(r1, r1size * 2 * sizeof(Int)) ;
                        Gnew = (double *) REALLOC(G2, r1size*4*sizeof(double)) ;
                        if (rnew != NULL)
                        {
                            r1 = rnew ;
                            G2 = Gnew ;
                            r1size *= 2 ;
                        }
                    }
                    G2[2*rrot_count] = 1 ;
                    G2[2*rrot_count+1] = 0 ;
                    r1[rrot_count++] = col ;
                }
            }

            fcol = ecol ;
            
            PRINT(("In Par row rot block\n")) ;
            /* Row rotations only R block */
            for (ccol = ecol+1 ; ccol < encol[srow] ; ccol++)
            {
                fcol = ccol ;
                cpk = cp[ccol] ;
                for (rindex = 0 ; rindex < rrot_count ; rindex++)
                {
                    c = G2[2*rindex] ;
                    s = G2[2*rindex+1] ;
                    row = r1[rindex] ;
                    PRINT(("Apply rowrot r-%ld-%ld c-%ld\n", row-1, row, ccol));
                    PRINT(("rot no %ld\n", rindex)) ;

                    da = A[cpk - (ccol-(row-1))] ;
                    db = A[cpk - (ccol-row)] ;
                    APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                    A[cpk - (ccol-(row-1))] = da ;
                    A[cpk - (ccol-row)] = db ;
                }
            }
            fcol++ ;

            PRINT(("In L block.\n")) ;

            rot_count = 0 ;
            if (fcol > n-1)
            {
                PRINT(("Stopping as row rots need not complete\n")) ;
                rrot_count = 0 ;
                continue ;
            }

            nexterow = ecol ;
            nextecol = -1 ;
            nextsrow = -1 ;
            PRINT(("rrot_count = %ld\n", rrot_count)) ;
            PRINT_INT_ARRAY (rrot_count, r1, "r1") ;

            for (rindex = 0 ; rindex < rrot_count ; rindex++)
            {
                c = G2[2*rindex] ;
                s = G2[2*rindex+1] ;
                row = r1[rindex] ;
                PRINT(("row = %ld\n", row)) ;
                ASSERT(row >= 0 && row < m, "row number out of bounds") ;
                lcol = encol[row] ;
                lcol1 = encol[row-1] ;

                PRINT(("rindex = %ld\n", rindex)) ;
                /*PRINT_INT_ARRAY (rrot_count, r1, "r1") ;*/
                PRINT(("lcol = %ld\n", lcol)) ;
                PRINT(("lcol1 = %ld\n", lcol1)) ;
                ASSERT(lcol >= 0 && lcol < n, "lcol out of bounds") ;
                ASSERT(lcol1 >= 0 && lcol1 < n, "lcol1 out of bounds") ;

                PRINT(("Apply rowrot r-%ld-%ld c-%ld-%ld\n", 
                            row-1, row, fcol, lcol)) ;
                PRINT(("rot no %ld\n", rindex)) ;

                /* Apply pending row rotations */
                for (cindex = fcol ; cindex <= lcol1 ; cindex++)
                {
                    cpk = cp[cindex] ;
                    ASSERT(cpk-(cindex-(row-1)) >= 0 && 
                        cpk-(cindex-(row-1)) < size, "A index out of bounds") ;
                    ASSERT(cpk-(cindex-(row)) >= 0 && 
                        cpk-(cindex-(row)) < size, "A index out of bounds") ;
                    da = A[cpk - (cindex-(row-1))] ;
                    db = A[cpk - (cindex-row)] ;
                    APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                    A[cpk - (cindex-(row-1))] = da ;
                    A[cpk - (cindex-row)] = db ;
                    fflush(stdout) ;
                }

                if (lcol1 < lcol)
                {
                    PRINT(("Creating fill \n"));
                    ASSERT(lcol1+1 == lcol, "More than one fill") ;
                    cpk = cp[lcol] ;
                    da = 0.0 ;
                    db = A[cpk - (lcol-row)] ;
                    PRINT(("cpk-(lcol-row) = %ld\n", cpk-(lcol-row)));
                    PRINT(("cpk=%ld lcol=%ld row=%ld\n", cpk, lcol, row));
                    ASSERT(cpk-(lcol-(row)) >= 0 && cpk-(lcol-(row)) < size, 
                        "A index out of ounds") ;
                    APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                    ws = da ;
                    A[cpk - (lcol-row)] = db ;

                    ASSERT(row !=0 , "row == 0 invalid here ") ;
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
                        PRINT(("chase corner %ld \n", lcol)) ;

                        if (ws != 0.0)
                        {
                            da = A[cp[lcol-1] - ((lcol-1)-(row-1))] ;
                            ASSERT(cp[lcol-1]-((lcol-1)-(row-1)) >= 0 && 
                                    cp[lcol-1]-((lcol-1)-(row-1)) < size, 
                                    "A index out of bounds") ;
                            db = ws ;
                            GIVENS(da, db, d, c, s) ;
                            A[cp[lcol-1] - ((lcol-1)-(row-1))] = d ; /* HERE */
                            ws = 0.0 ;
                            PRINT(("F&A Givens r-%ld-%ld c-%ld-%ld\n", row-1, 
                                        ecol,
                            lcol-1, lcol)) ;

                            cpk = cp[lcol] ;
                            cpk1 = cp[lcol-1] ;
                            for (rowindex = row ; rowindex <= ecol ; rowindex++)
                            {
                                da = A[cpk1 - ((lcol-1)-rowindex)] ;
                                db = A[cpk - (lcol-rowindex)] ;
                                APPLY_GIVENS(da, db, c, s, dtemp, conjs) ;
                                A[cpk1 - ((lcol-1)-rowindex)] = da ;
                                A[cpk - (lcol-rowindex)] = db ;
                            }

                            if (sfact->V != NULL)
                            {
                                APPLY_GIVENS_TO_COL(sfact->V, n, lcol-1, lcol,
                                    0, n-1, c, s, tmpi, da, db, dtemp, n) ;
                            }


                            if (rot_count == c1size)
                            {
                            PRINT(("Realloc of c1 and G1 ***************")) ;
                                ASSERT(opt == 2, "opt != 2") ;
                                cnew = (Int *) REALLOC 
                                        (c1, c1size * 2 * sizeof(Int)) ;
                                Gnew = (double *) REALLOC
                                        (G1, c1size*4*sizeof(double)) ;
                                if (cnew != NULL && Gnew != NULL)
                                {
                                    c1 = cnew ;
                                    G1 = Gnew ;
                                    c1size *= 2 ;
                                }
                                else
                                {
                                    ASSERT(0, "Out of Memory") ;
                                }
                            }
                            ASSERT((2 * rot_count+1) < c1size*2  && 
                                    (2 * rot_count+1) > 0, "invalid rot_count");
                            G1[2*rot_count] = c ;
                            G1[2*rot_count+1] = s ;
                            c1[rot_count++] = lcol ;
                        }
                        else
                        {
                            if (rot_count == c1size)
                            {
                                PRINT(("Realloc of c1 and G1 ************")) ;
                                ASSERT(opt == 2, "opt != 2") ;
                                cnew = (Int *) REALLOC 
                                        (c1, c1size * 2 * sizeof(Int)) ;
                                Gnew = (double *) REALLOC
                                        (G1, c1size*4*sizeof(double)) ;
                                if (cnew != NULL && Gnew != NULL)
                                {
                                    c1 = cnew ;
                                    G1 = Gnew ;
                                    c1size *= 2 ;
                                }
                                else
                                {
                                    ASSERT(0, "Out of Memory") ;
                                }
                            }
                            ASSERT((2 * rot_count+1) < c1size*2  && 
                            (2 * rot_count+1) > 0, "invalid rot_count ") ;
                            G1[2*rot_count] = 1 ;
                            G1[2*rot_count+1] = 0 ;
                            c1[rot_count++] = lcol ;
                        }

                        /* set trow and end encol to old values */
                        encol[row-1] = lcol1 ;
                        trow[lcol] = row ;
                        PRINT(("lcol= "ID" trow[lcol]="ID" bw = "ID"\n",
                            lcol, trow[lcol], bw)) ;
                        ASSERT(lcol - trow[lcol] >= bw, "trow[lcol] wrong") ;
                        /* The fill is chasable, min_trow will remain the same*/

                    }
                    else
                    {

                        PRINT(("row-1= "ID" encol[row-1]="ID" bw = "ID"\n",
                            row-1, encol[row-1], bw)) ;
                        if (encol[row-1] != n-1)
                        {
                            ASSERT(encol[row-1]-(row-1) >= bw ,
                                "encol[k] wrong") ;
                        }
                        PRINT(("lcol= "ID" trow[lcol]="ID" bw = "ID"\n", lcol,
                            trow[lcol], bw)) ;
                        ASSERT(lcol - trow[lcol] >= bw, "trow[lcol] wrong") ;

                        /* Fill that cannot be chased */
                        PRINT(("%ld is not a corner\n", lcol)) ;

                        /* TBD : Need to fix linked list stuff */
                        if ((cp[lcol] - (lcol-(row-1))) > cp[prev[lcol]])
                        /*if ((cp[lcol] - (lcol-(row-1))) > cp[lcol-1])*/
                        {
                            /* There is enough space */
                            A[cp[lcol] - (lcol-(row-1))] = ws ;
                        }
                        else
                        {
                            PRINT(("Before the move \n")) ;
                            PRINT_SKYLINE (sym, sky) ;
                            /* Need to move column lcol. Additional space is in
                             * A[cp[tail]..cp[n+1]].
                             * */

#if 0
                            PRINT(("trow[%d] = %d\n", lcol, trow[lcol])) ;
                            PRINT(("min_trow[%d] = %d\n", lcol, 
                                        min_trow[lcol]));
                            PRINT(("lcol = %d\n", lcol)) ;
                            PRINT(("cp[%d] = %d\n", lcol-1, cp[lcol-1])) ;
                            PRINT(("cp[%d] = %d\n", lcol, cp[lcol])) ;
                            PRINT(("row = %d\n", row)) ;
#endif

                            ASSERT (opt == 3, "opt != 3") ;

                            avail = cp[n] - cp[tail] ;
                            lcol_size = lcol - trow[lcol] + 1 ;
                            /* TBD : Change 5 based on max column ht later.
                             * No need to add 1 for ws, as trow[lcol] is 
                             * changed already. */
                            required = lcol_size + 5 ;
                            PRINT(("required = %ld\n", required)) ;

                            PRINT(("cp[%ld]=%ld, prev[%ld]=%ld cp[%ld]=%ld\n", 
                            lcol, cp[lcol], lcol, prev[lcol], prev[lcol], 
                            cp[prev[lcol]])) ;

                            PRINT(("new required = %ld, avail=%ld\n", 
                                required, avail)) ;
                            PRINT(("cp[tail]=%ld\n", cp[tail])) ;
                            ASSERT(required > 0, "reqd <= 0") ;

                            if (avail >= required)
                            {
                                ASSERT((cp[tail]+6) >= 0 && (cp[tail]+6) < size,
                                "A index out of ounds") ;

                                A[cp[tail] + 6] = ws ; 
                                tindex = cp[tail] + 7 ;
                                PRINT(("size = %ld cp[n]=%ld \n", size, cp[n])) ;
                                ASSERT(size == cp[n]+1, "size != cp[n]+1") ;
                                /* We add 2 instead of 1 bcos lcol_size already
                                 * includes ws as trow changed 
                                 * */
                                cindex = cp[lcol] - lcol_size + 2 ;
                                for ( ; cindex <= cp[lcol] ; cindex++)
                                {
                                    PRINT(("tindex =%ld\n", tindex)) ;
                                    ASSERT(tindex >= 0 && tindex < size, 
                                        "A index out of bounds") ;
                                    ASSERT(cindex >= 0 && cindex < size, 
                                        "A index out of bounds") ;

                                    A[tindex++] = A[cindex] ;
                                }

                                ASSERT(lcol >= 0 && lcol < n, 
                                        "cp index out of ounds") ;
                                cp[lcol] = tindex - 1 ; /* HERE */

                                PRINT(("Moving column %ld \n", lcol)) ;
                                PRINT(("head = %ld, tail = %ld\n", head, tail));
                                PRINT(("lcol prev=%ld, next=%ld\n", 
                                        prev[lcol], next[lcol])) ;
                                PRINT(("ctail - cnext = %ld\n", next[tail])) ;

                                delete_column(lcol, next, prev, &head, &tail) ;
                                add_column(lcol, next, prev, &head, &tail) ;

                                PRINT(("After change \n")) ;
                                PRINT(("head = %ld, tail = %ld\n", head, tail));
                                PRINT(("lcol prev=%ld, next=%ld\n", 
                                        prev[lcol], next[lcol])) ;
                                PRINT(("ctail - cnext = %ld\n", next[tail])) ;
                                ASSERT(tail == lcol, "tail != lcol") ;
                            }
                            else
                            {
                                /* Not enough space create more space and then
                                 * move the column. 
                                 * TBD : Should this 3 * n be elbow too ?
                                 * */
                                Gnew = (double *) REALLOC(A, ((cp[n]+1)+3*n)*
                                                            sizeof(double)) ;
                                if (Gnew != NULL)
                                {
                                    A = Gnew ;
                                    sky->A = Gnew ;
                                }
                                else
                                {
                                    PRINT(("Out of memory \n")) ;
                                    /* TBD :Should call pack here */
                                }
                                for (tindex = size ; tindex < size + 3 *n ;
                                        tindex++)
                                {
                                    A[tindex] = 0.0 ;
                                }
                                realloc_cnt++ ;
                                size = size + 3 * n ;
                                /* Actual index not the size , cp[tail] != cp[n]
                                 * */
                                cp[n] = cp[n] + 3 * n ; /* HERE */

                                ASSERT((cp[tail]+6) >= 0 && (cp[tail]+6) < size,
                                "A index out of ounds") ;

                                A[cp[tail]+6] = ws ;
                                tindex = cp[tail] + 7 ;

                                /* We add 2 instead of 1 bcos lcol_size already
                                 * includes ws as trow changed 
                                 * */
                                cindex = cp[lcol] - lcol_size + 2 ;
                                for ( ; cindex <= cp[lcol] ; cindex++)
                                {
                                    ASSERT(tindex >= 0 && tindex < size, 
                                        "A index out of ounds") ;
                                    ASSERT(cindex >= 0 && cindex < size, 
                                        "A index out of ounds") ;
                                    A[tindex++] = A[cindex] ;
                                }

                                ASSERT(lcol >= 0 && lcol < n, 
                                        "cp index out of ounds") ;
                                cp[lcol] = tindex - 1 ;

                                delete_column(lcol, next, prev, &head, &tail) ;
                                add_column(lcol, next, prev, &head, &tail) ;

                                ASSERT(tail == lcol, "tail != lcol") ;

                            }
                            PRINT(("After the move \n")) ;
                            PRINT_SKYLINE (sym, sky) ;
                        }

                        /*min_trow[lcol] = MIN(min_trow[lcol], trow[lcol]) ;*/
                    }
                }
                else
                {
                    PRINT(("%ld has no fill\n", lcol)) ;
                    fflush(stdout) ;
                }
            }
            rrot_count = 0 ;
            srow = nextsrow ;
            erow = nexterow ;
            ecol = nextecol ;
            PRINT(("nextsrow = %ld\n", srow)) ;
            PRINT(("nexterow = %ld\n", erow)) ;
            PRINT(("nextecol = %ld\n", ecol)) ;
            /*PRINT_SKYLINE (sym, sky) ;*/
        }
        PRINT(("Completed the chase\n")) ;

        }
        PRINT(("top3=%ld\n", top3)) ;
        /* pop pending corners to rmc stack list */
        while (nlist3 > top3)
        {
            PRINT(("Adding %ld to right stack\n", list3[nlist3-1])) ;
            list[nlist++] = list3[nlist3-1] ;
            ASSERT (list3[nlist3-1] >= 0, "list3[nlist3-1] < 0") ;
            list3[--nlist3] = -1 ;
        }
        nlist3 = 0 ;
        top3 = 0 ;
        PRINT_SKYLINE_SVD(sym, sky, sgood) ;

    }
    scom->realloc_cnt = realloc_cnt ; 
    /*mexPrintf("Realloc count = %ld \n", realloc_cnt) ;*/

    PRINT_SKYLINE_SVD(sym, sky, sgood) ;

    /* ---------------------------------------------------------------------- */
    /* TODO describe me */
    /* ---------------------------------------------------------------------- */

    if (scom->bw_est != 1)
    {
        /* copy the matrix into a band and call the band code explicitly */

        PRINT_SKYLINE_SVD(sym, sky, sgood) ;

        band_A = (double *) CALLOC ((bw+1) * n , sizeof(double)) ;
        ldab = bw + 1 ;
        obu = bw ;

        for (k = 0 ; k < n ; k++)
        {
            /*PRINT(("Column %ld", k)) ;*/
            for (i = trow[k] ; i <= k ; i++)
            {
                /*PRINT((" %0.4f", A[cp[k] - (k-i)])) ;*/
                band_A[INDEX(i, k)] = A[cp[k] - (k-i)] ;
            }
            /*PRINT(("\n")) ;
            fflush(stdout) ;*/
        }

        sfact->b2[0] = 0.0 ;
        piro_band_reduce_drl(NULL, n, n, 0, 0, bw, band_A, ldab, sfact->b1, 
                        sfact->b2+1, NULL, 0, NULL, 0, NULL, 0, NULL, 0) ;


        FREE(band_A) ;

    }
    else
    {

        /* Need to assign bidiagonal matrix to return : TBD */
        sfact->b1[0] = A[cp[0]] ;
        sfact->b2[0] = 0.0 ;
        for (i = 1 ; i < MIN(m, n) ; i++)
        {
            sfact->b1[i] = A[cp[i]] ;
            sfact->b2[i] = A[cp[i]-1] ;
        }
    }

        /*PRINT(("nlist = %ld \n", nlist)) ;
        PRINT_INT_ARRAY (nlist, list, "list") ;
        PRINT_INT_ARRAY (n+1, trow, "trow") ;
        PRINT_INT_ARRAY (n, encol, "encol") ;
        PRINT_INT_ARRAY (n, min_trow, "min_trow") ;*/

    FREE(intmem) ;
    FREE(tempcol) ;
    FREE(sgood) ;
    if (opt == 2)
    {
        FREE(r1) ;
        FREE(c1) ;
        FREE(G1) ;
        FREE(G2) ;
    }

    PRINT(("Leaving sky_reduce \n")) ;
    return (0) ;
}
