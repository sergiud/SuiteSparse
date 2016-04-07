#include <stdlib.h>
#include "genband_util.h"
#include "genband_internal.h"
#include "blksky.h"

/* ========================================================================== */

Int sky_set_common              /* TODO: why return Int? */
(
    sky_common *scom
)
{
    double safemin, eps ;
    Int esfmn2 ;


    /* TBD : Should adjust for float */
    safemin = GENBAND_DBL_MIN ;
    eps = GENBAND_EPS/2 ;
    esfmn2  = log(safemin/eps) / log(2) / 2 ;
    scom->safemin2 = pow(2, esfmn2) ; 
    scom->safemx2 = 1/(scom->safemin2) ;

    PRINT (("sky_set_common:\n")) ;
    PRINT (("safemin  %g\n", safemin)) ;
    PRINT (("eps      %g\n", eps)) ;
    PRINT (("esfmn2   %ld\n", esfmn2)) ;
    PRINT (("safemin2 %g\n", scom->safemin2)) ;
    PRINT (("safemx2  %g\n", scom->safemx2)) ;

    scom->opt = 0 ;
    scom->cspace = 0 ;
    scom->elbow = 0 ;
    scom->rcount = 0 ;
    scom->realloc_cnt = 0 ;
    scom->iwork = NULL ;
    scom->dwork = NULL ;
    scom->bw_est = 0 ;

    return (0) ;
}

/* ========================================================================== */

Int sky_allocate_symbolic       /* TODO why return Int? */
(
    sky_common *scom, 
    sky_symbolic *sym 
)
{
    Int n ;

    /* TODO: n should be an input parameter, not in sym->n already */
    n = sym->n ;


    /* NOTE: sym->cp does not exist */

    /* TODO: are there parts of sym not initialized? */

    /* TODO if n is zero, many terms below are NULL.  This is bad */

    sym->trow = (Int *) MALLOC((n+1) * sizeof(Int)) ;
    sym->encol = (Int *) MALLOC(n * sizeof(Int)) ;
    sym->prev = (Int *) MALLOC(n * sizeof(Int)) ;
    sym->next = (Int *) MALLOC(n * sizeof(Int)) ;
    if (scom->opt != 3)
    {
        /* TODO Help!  min_trow is not described anywhere */
        sym->min_trow = (Int *) MALLOC(n * sizeof(Int)) ;
    }
    else
    {
        sym->min_trow = NULL ;
    }

    /* TODO what about giv_cnt, swap_cnt, flp_cnt, etc? */

    return (0) ;
}

/* ========================================================================== */

Int sky_free_symbolic           /* TODO why return Int? */
(
    sky_symbolic *sym 
)
{
    FREE(sym->trow) ;
    FREE(sym->encol) ;
    FREE(sym->prev) ;
    FREE(sym->next) ;
    if (sym->min_trow != NULL) FREE(sym->min_trow) ;
    return (0) ;
}

/* ========================================================================== */


Int sky_allocate_skyline        /* TODO why return Int? */
(
    sky_common *scom,
    sky_skyline *sky,
    sky_symbolic *sym
)
{
    Int clen, maxht ;
    Int cspace, elbow ;
    Int csize ;
    Int *m_trow ;
    Int i ;
    Int n ;

    /* TODO: what is m_trow and why is it different than trow ? */

    m_trow = (sym->min_trow != NULL) ? sym->min_trow : sym->trow ;
    cspace = scom->cspace ;
    elbow = scom->elbow ;

    n = sky->n ;
    csize = 0 ;

    maxht = 0 ;
    /* Allocate column pointers in the skyline */
    sky->cp = (Int *) MALLOC((n+1) * sizeof(Int)) ;

    for (i = 0 ; i < n ; i++)
    {
        ASSERT(m_trow[i] >= 0, "m_trow[i] < 0") ;
        clen = i - m_trow[i] + 1 ;

        /*
        printf ("i %g m_trow %g (%g %g)clen %g\n",
            (double) i, (double) m_trow [i],
            (double) ((sym->min_trow == NULL) ? (-999) : sym->min_trow [i]),
            (double) (sym->trow [i]),
            (double) clen) ;
        */
            
#if 0
        maxht = MAX(maxht, clen) ; /* Check this for opt == 1 TBD !!!!! */
        clen = MIN(clen+cspace, maxht) ;
#endif
        clen = clen + cspace ;
        PRINT(("clen %ld \n", clen)) ;

        csize = csize + clen ;
        PRINT(("i=%ld csize=%ld \n", i, csize)) ;
        sky->cp[i] = csize - 1 ;     /* The index, not the actual size */
        PRINT(("csize "ID"\n", csize)) ;
        ASSERT (csize > 0, "csize <= 0") ;
    }

    sky->cp[n] = csize + elbow - 1 ; /* The index not the actual size */

    ASSERT(csize > 0, "csize <= 0") ;
    PRINT(("Allocating "ID" bytes\n", csize+elbow)) ;

    /* TODO: Ugh.  Do not use CALLOC */
    sky->A = (double *) CALLOC (csize+elbow , sizeof(double)) ;
    return (0) ;
}

/* ========================================================================== */

Int sky_free_skyline                /* TODO why return Int? */
(
    sky_skyline *sky
)
{
    FREE (sky->cp) ;    /* TODO: does this set sky->cp to NULL as it should? */
    FREE (sky->A) ; 
    return (0) ;
}

/* ========================================================================== */

/* Need sky->cp to be set correctly */
Int sky_sparse_to_skyline           /* TODO why return Int? */
(
    sky_skyline *sky,   /* TODO comment me */
    Int *Ap,
    Int *Ai,
    double *Ax
)
{
    Int p, i, k ;
    Int *cp ;
    double *A ;

    cp = sky->cp ;
    A = sky->A ;

    /* note that this does not require A to have sorted row indices */

    /* TODO: modify this so that A need not be zero on input (see CALLOC
       above) */

    /* TODO: what about duplicate entries?  Either detect them and flag an
       error, or sum them up */

    for (k = 0 ; k < sky->n ; k++)
    {
        for (p = Ap[k] ; p < Ap[k+1] ; p++)
        {
            i = Ai [p] ;
            A [cp [k] - (k - i)] = Ax [p] ; 
        }
    }
    return (0) ;
}
