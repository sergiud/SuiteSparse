/* Requires CXSparse */

#include "cs.h"

#include "genband_util.h"
#include "genband_internal.h"
#include "blksky.h"

static void blksky_reduce_all
(
    cs_dl *R 
)
{
    Int m, n ;
    double *Ax ;
    Int *Ap, *Ai ;
    Int csize ;
    Int nnz ;
    Int i ;
    sky_common sc,  *scom ;
    sky_symbolic sm,  *ssym ;
    sky_skyline ss, *ssky ;
    sky_factor sf, *sfact ;

    scom = &sc ;
    sky_set_common (scom) ;

    ssym = &sm ;
    ssky = &ss ;
    sfact = &sf ;

    m = R->m ;
    n = R->n ;

    ssym->n = n ;

    ssky->m = m ;
    ssky->n = n ;

    sfact = &sf ;
    sfact->b1 = (double *) malloc(MIN(m, n) * sizeof(double)) ;
    sfact->b2 = (double *) malloc(MIN(m, n) * sizeof(double)) ;
    sfact->U = NULL ;
    sfact->V = NULL ;

    Ap = (Int *) R->p ;
    Ai = (Int *) R->i ;
    Ax = R->x ;
    nnz = Ap[n] ;

#if 0
    scom->opt = 3 ; /* TBD : for now */

    /*if (nargin < 2)
    {
        scom->opt = 1 ;
    }
    else
    {
        scom->opt = (Int) mxGetScalar(prhs[1]) ;
    }

    if (scom->opt == 3)
    {*/
        scom->cspace = 5 ;
        scom->elbow = 3 * n ;
    /*}
    else
    {
        scom->cspace = 0 ;
        scom->elbow = 0 ;
    }*/
#endif

    scom->opt = 1 ;
    scom->cspace = 0 ;
    scom->elbow = 0 ;

    scom->bw_est = 1 ;

    sky_allocate_symbolic(scom, ssym) ;

    blksky_find_symbolic(scom, Ap, Ai, Ax, m, n, ssym) ;

    sky_allocate_skyline(scom, ssky, ssym) ;

    csize = ssky->cp[n] ;
    /*if (nnz != 0)
    {
        printf("nnz=%d, csize=%d, csize/nnz = %d\n", nnz, csize, csize/nnz) ;
    }
    else
    {
        printf("nnz=%d, csize=%d, \n", nnz, csize) ;
    }*/

    sky_sparse_to_skyline(ssky, Ap, Ai, Ax) ;

    if (scom->opt != 2)
    {
        scom->rcount = 32 ;
        scom->dwork = (double *) malloc( 4 * scom->rcount * sizeof(double)) ;
        scom->iwork = (Int *) malloc( 2 * scom->rcount * sizeof(Int)) ;
    }
    else
    {
        scom->rcount = -1 ;
        scom->dwork = NULL ;
        scom->iwork = NULL ;
    }
    sky_reduce (scom, ssym, ssky, sfact) ;

    /*plhs[0] = mxCreateDoubleMatrix(MIN(m, n), 1, mxREAL) ;
    plhs[1] = mxCreateDoubleMatrix(MIN(m, n), 1, mxREAL) ;
    rb1 = mxGetPr(plhs[0]) ;
    rb2 = mxGetPr(plhs[1]) ;*/
    /* copy b1 and b2 to rb1 and rb2 */
    for (i = 0 ; i < MIN(m, n) ; i++)
    {
        /*rb1[i] = sfact->b1[i] ;
        rb2[i] = sfact->b2[i] ;*/
        printf("b1["ID"] = %g, b2["ID"] = %g\n",
            i, sfact->b1[i], i, sfact->b2[i]) ;
    }

    free(sfact->b1) ;
    free(sfact->b2) ;
    if (scom->opt != 2)
    {
        free(scom->dwork) ;
        free(scom->iwork) ;
    }

    sky_free_skyline (ssky) ; 
    sky_free_symbolic (ssym) ;

}

int main(void)
{
    cs_dl *T, *R ;
    FILE *fid ;
    Int m, n, i ;

    printf ("sizeof Int %g\n", (double) sizeof (Int)) ;

    /*
    fid = fopen("Rmat", "r") ;
    if (fid == NULL) return ;
    */

    /* get the matrix from stdin */
    fid = stdin ;

    /* first line has # of rows and columns */
    m = 0 ;
    n = 0 ;
    i = fscanf (fid, ID" "ID"\n", &m, &n) ;
    if (i < 2) printf ("warning: file may be corrupted\n") ;

    /* load the rest of the triplets */
    T = cs_dl_load(fid) ;
    T->m = MAX (T->m, m) ;
    T->n = MAX (T->n, n) ;
    printf ("T is "ID" by "ID"\n", T->m, T->n) ;

    /* convert the triplets into a compressed-column matrix */
    R = cs_dl_compress (T) ;

    /* free the triplet matrix */
    T = cs_dl_spfree (T) ;

    /* print the first part of R */
    cs_dl_print (R, 1) ;

    /* reduce R to bidiagonal form */
    blksky_reduce_all (R) ;

    /* free R */
    R = cs_dl_spfree (R) ;

    return 1 ;
}
