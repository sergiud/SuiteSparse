//------------------------------------------------------------------------------
// CHOLMOD/Supernodal/t_cholmod_super_solve: template for cholmod_super_solve
//------------------------------------------------------------------------------

// CHOLMOD/Supernodal Module.  Copyright (C) 2005-2022, Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

/* Template routine for cholmod_super_solve.  Supports real or complex L. */

#include "cholmod_template.h"

static void TEMPLATE (cholmod_super_lsolve)
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to use for the forward solve */
    /* ---- output ---- */
    cholmod_dense *X,	/* b on input, solution to Lx=b on output */
    /* ---- workspace ---- */
    cholmod_dense *E,	/* workspace of size nrhs*(L->maxesize) */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Lx, *Xx, *Ex ;
    double minus_one [2], one [2] ;
    Int *Lpi, *Lpx, *Ls, *Super ;
    Int nsuper, k1, k2, psi, psend, psx, nsrow, nscol, ii, s,
	nsrow2, n, ps2, j, i, d, nrhs ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrhs = X->ncol ;
    Ex = E->x ;
    Xx = X->x ;
    n = L->n ;
    d = X->d ;

    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;
    Lx = L->x ;
    minus_one [0] = -1.0 ;
    minus_one [1] = 0 ;
    one [0] = 1.0 ;
    one [1] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* solve Lx=b */
    /* ---------------------------------------------------------------------- */

    if (nrhs == 1)
    {

	for (s = 0 ; s < nsuper ; s++)
	{
	    k1 = Super [s] ;
	    k2 = Super [s+1] ;
	    psi = Lpi [s] ;
	    psend = Lpi [s+1] ;
	    psx = Lpx [s] ;
	    nsrow = psend - psi ;
	    nscol = k2 - k1 ;
	    nsrow2 = nsrow - nscol ;
	    ps2 = psi + nscol ;
	    ASSERT ((size_t) nsrow2 <= L->maxesize) ;

	    /* L1 is nscol-by-nscol, lower triangular with non-unit diagonal.
	     * L2 is nsrow2-by-nscol.  L1 and L2 have leading dimension of
	     * nsrow.  x1 is nscol-by-nsrow, with leading dimension n.
	     * E is nsrow2-by-1, with leading dimension nsrow2.
	     */

	    /* gather X into E */
	    for (ii = 0 ; ii < nsrow2 ; ii++)
	    {
		/* Ex [ii] = Xx [Ls [ps2 + ii]] ; */
		ASSIGN (Ex,-,ii, Xx,-,Ls [ps2 + ii]) ;
	    }

#ifdef REAL

	    /* solve L1*x1 (that is, x1 = L1\x1) */
	    SUITESPARSE_BLAS_dtrsv ("L", "N", "N",
		nscol,			    /* N:       L1 is nscol-by-nscol */
		Lx + ENTRY_SIZE*psx, nsrow, /* A, LDA:  L1 */
		Xx + ENTRY_SIZE*k1, 1,      /* X, INCX: x1 */
                Common->blas_ok) ;

	    /* E = E - L2*x1 */
	    SUITESPARSE_BLAS_dgemv ("N",
		nsrow2, nscol,		    /* M, N:    L2 is nsrow2-by-nscol */
		minus_one,		    /* ALPHA:   -1 */
		Lx + ENTRY_SIZE*(psx + nscol),   /* A, LDA:  L2 */
		nsrow,
		Xx + ENTRY_SIZE*k1, 1,	    /* X, INCX: x1 */
		one,			    /* BETA:    1 */
		Ex, 1,                      /* Y, INCY: E */
                Common->blas_ok) ;

#else

	    /* solve L1*x1 (that is, x1 = L1\x1) */
	    SUITESPARSE_BLAS_ztrsv ("L", "N", "N",
		nscol,			    /* N:       L1 is nscol-by-nscol */
		Lx + ENTRY_SIZE*psx, nsrow, /* A, LDA:  L1 */
		Xx + ENTRY_SIZE*k1, 1,      /* X, INCX: x1 */
                Common->blas_ok) ;

	    /* E = E - L2*x1 */
	    SUITESPARSE_BLAS_zgemv ("N",
		nsrow2, nscol,		    /* M, N:    L2 is nsrow2-by-nscol */
		minus_one,		    /* ALPHA:   -1 */
		Lx + ENTRY_SIZE*(psx + nscol),   /* A, LDA:  L2 */
		nsrow,
		Xx + ENTRY_SIZE*k1, 1,	    /* X, INCX: x1 */
		one,			    /* BETA:    1 */
		Ex, 1,                      /* Y, INCY: E */
                Common->blas_ok) ;

#endif

	    /* scatter E back into X */
	    for (ii = 0 ; ii < nsrow2 ; ii++)
	    {
		/* Xx [Ls [ps2 + ii]] = Ex [ii] ; */
		ASSIGN (Xx,-,Ls [ps2 + ii], Ex,-,ii) ;
	    }
	}
    }
    else
    {

	for (s = 0 ; s < nsuper ; s++)
	{
	    k1 = Super [s] ;
	    k2 = Super [s+1] ;
	    psi = Lpi [s] ;
	    psend = Lpi [s+1] ;
	    psx = Lpx [s] ;
	    nsrow = psend - psi ;
	    nscol = k2 - k1 ;
	    nsrow2 = nsrow - nscol ;
	    ps2 = psi + nscol ;
	    ASSERT ((size_t) nsrow2 <= L->maxesize) ;

	    /* E is nsrow2-by-nrhs, with leading dimension nsrow2. */

	    /* gather X into E */
	    for (ii = 0 ; ii < nsrow2 ; ii++)
	    {
		i = Ls [ps2 + ii] ;
		for (j = 0 ; j < nrhs ; j++)
		{
		    /* Ex [ii + j*nsrow2] = Xx [i + j*d] ; */
		    ASSIGN (Ex,-,ii+j*nsrow2, Xx,-,i+j*d) ;
		}
	    }

#ifdef REAL

	    /* solve L1*x1 */
	    SUITESPARSE_BLAS_dtrsm ("L", "L", "N", "N",
		nscol, nrhs,			/* M, N: x1 is nscol-by-nrhs */
		one,				/* ALPHA:  1 */
		Lx + ENTRY_SIZE*psx, nsrow,	/* A, LDA: L1 */
		Xx + ENTRY_SIZE*k1, d,          /* B, LDB: x1 */
                Common->blas_ok) ;

	    /* E = E - L2*x1 */
	    if (nsrow2 > 0)
	    {
		SUITESPARSE_BLAS_dgemm ("N", "N",
		    nsrow2, nrhs, nscol,	    /* M, N, K */
		    minus_one,			    /* ALPHA:  -1 */
		    Lx + ENTRY_SIZE*(psx + nscol),  /* A, LDA: L2 */
		    nsrow,
		    Xx + ENTRY_SIZE*k1, d,	    /* B, LDB: X1 */
		    one,			    /* BETA:   1 */
		    Ex, nsrow2,                     /* C, LDC: E */
                    Common->blas_ok) ;
	    }

#else

	    /* solve L1*x1 */
	    SUITESPARSE_BLAS_ztrsm ("L", "L", "N", "N",
		nscol, nrhs,			/* M, N: x1 is nscol-by-nrhs */
		one,				/* ALPHA:  1 */
		Lx + ENTRY_SIZE*psx, nsrow,	/* A, LDA: L1 */
		Xx + ENTRY_SIZE*k1, d,          /* B, LDB: x1 */
                Common->blas_ok) ;

	    /* E = E - L2*x1 */
	    if (nsrow2 > 0)
	    {
		SUITESPARSE_BLAS_zgemm ("N", "N",
		    nsrow2, nrhs, nscol,	    /* M, N, K */
		    minus_one,			    /* ALPHA:  -1 */
		    Lx + ENTRY_SIZE*(psx + nscol),  /* A, LDA: L2 */
		    nsrow,
		    Xx + ENTRY_SIZE*k1, d,	    /* B, LDB: X1 */
		    one,			    /* BETA:   1 */
		    Ex, nsrow2,                     /* C, LDC: E */
                    Common->blas_ok) ;
	    }

#endif

	    /* scatter E back into X */
	    for (ii = 0 ; ii < nsrow2 ; ii++)
	    {
		i = Ls [ps2 + ii] ;
		for (j = 0 ; j < nrhs ; j++)
		{
		    /* Xx [i + j*d] = Ex [ii + j*nsrow2] ; */
		    ASSIGN (Xx,-,i+j*d, Ex,-,ii+j*nsrow2) ;
		}
	    }
	}
    }
}


static void TEMPLATE (cholmod_super_ltsolve)
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to use for the forward solve */
    /* ---- output ---- */
    cholmod_dense *X,	/* b on input, solution to Lx=b on output */
    /* ---- workspace ---- */
    cholmod_dense *E,	/* workspace of size nrhs*(L->maxesize) */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Lx, *Xx, *Ex ;
    double minus_one [2], one [2] ;
    Int *Lpi, *Lpx, *Ls, *Super ;
    Int nsuper, k1, k2, psi, psend, psx, nsrow, nscol, ii, s,
	nsrow2, n, ps2, j, i, d, nrhs ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrhs = X->ncol ;
    Ex = E->x ;
    Xx = X->x ;
    n = L->n ;
    d = X->d ;

    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;
    Lx = L->x ;
    minus_one [0] = -1.0 ;
    minus_one [1] = 0 ;
    one [0] = 1.0 ;
    one [1] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* solve L'x=b */
    /* ---------------------------------------------------------------------- */

    if (nrhs == 1)
    {

	for (s = nsuper-1 ; s >= 0 ; s--)
	{
	    k1 = Super [s] ;
	    k2 = Super [s+1] ;
	    psi = Lpi [s] ;
	    psend = Lpi [s+1] ;
	    psx = Lpx [s] ;
	    nsrow = psend - psi ;
	    nscol = k2 - k1 ;
	    nsrow2 = nsrow - nscol ;
	    ps2 = psi + nscol ;
	    ASSERT ((size_t) nsrow2 <= L->maxesize) ;

	    /* L1 is nscol-by-nscol, lower triangular with non-unit diagonal.
	     * L2 is nsrow2-by-nscol.  L1 and L2 have leading dimension of
	     * nsrow.  x1 is nscol-by-nsrow, with leading dimension n.
	     * E is nsrow2-by-1, with leading dimension nsrow2.
	     */

	    /* gather X into E */
	    for (ii = 0 ; ii < nsrow2 ; ii++)
	    {
		/* Ex [ii] = Xx [Ls [ps2 + ii]] ; */
		ASSIGN (Ex,-,ii, Xx,-,Ls [ps2 + ii]) ;
	    }

#ifdef REAL

	    /* x1 = x1 - L2'*E */
	    SUITESPARSE_BLAS_dgemv ("C",
		nsrow2, nscol,		    /* M, N: L2 is nsrow2-by-nscol */
		minus_one,		    /* ALPHA:   -1 */
		Lx + ENTRY_SIZE*(psx + nscol),   /* A, LDA:  L2 */
		nsrow,
		Ex, 1,			    /* X, INCX: Ex */
		one,			    /* BETA:    1 */
		Xx + ENTRY_SIZE*k1, 1,      /* Y, INCY: x1 */
                Common->blas_ok) ;

	    /* solve L1'*x1 */
	    SUITESPARSE_BLAS_dtrsv ("L", "C", "N",
		nscol,			    /* N:	L1 is nscol-by-nscol */
		Lx + ENTRY_SIZE*psx, nsrow,	    /* A, LDA:  L1 */
		Xx + ENTRY_SIZE*k1, 1,              /* X, INCX: x1 */
                Common->blas_ok) ;

#else

	    /* x1 = x1 - L2'*E */
	    SUITESPARSE_BLAS_zgemv ("C",
		nsrow2, nscol,		    /* M, N: L2 is nsrow2-by-nscol */
		minus_one,		    /* ALPHA:   -1 */
		Lx + ENTRY_SIZE*(psx + nscol),   /* A, LDA:  L2 */
		nsrow,
		Ex, 1,			    /* X, INCX: Ex */
		one,			    /* BETA:    1 */
		Xx + ENTRY_SIZE*k1, 1,      /* Y, INCY: x1 */
                Common->blas_ok) ;

	    /* solve L1'*x1 */
	    SUITESPARSE_BLAS_ztrsv ("L", "C", "N",
		nscol,			    /* N:	L1 is nscol-by-nscol */
		Lx + ENTRY_SIZE*psx, nsrow,	    /* A, LDA:  L1 */
		Xx + ENTRY_SIZE*k1, 1,              /* X, INCX: x1 */
                Common->blas_ok) ;

#endif

	}
    }
    else
    {

	for (s = nsuper-1 ; s >= 0 ; s--)
	{
	    k1 = Super [s] ;
	    k2 = Super [s+1] ;
	    psi = Lpi [s] ;
	    psend = Lpi [s+1] ;
	    psx = Lpx [s] ;
	    nsrow = psend - psi ;
	    nscol = k2 - k1 ;
	    nsrow2 = nsrow - nscol ;
	    ps2 = psi + nscol ;
	    ASSERT ((size_t) nsrow2 <= L->maxesize) ;

	    /* E is nsrow2-by-nrhs, with leading dimension nsrow2. */

	    /* gather X into E */
	    for (ii = 0 ; ii < nsrow2 ; ii++)
	    {
		i = Ls [ps2 + ii] ;
		for (j = 0 ; j < nrhs ; j++)
		{
		    /* Ex [ii + j*nsrow2] = Xx [i + j*d] ; */
		    ASSIGN (Ex,-,ii+j*nsrow2, Xx,-,i+j*d) ;
		}
	    }

#ifdef REAL

	    /* x1 = x1 - L2'*E */
	    if (nsrow2 > 0)
	    {
		SUITESPARSE_BLAS_dgemm ("C", "N",
		    nscol, nrhs, nsrow2,	/* M, N, K */
		    minus_one,			/* ALPHA:  -1 */
		    Lx + ENTRY_SIZE*(psx + nscol),  /* A, LDA: L2 */
		    nsrow,
		    Ex, nsrow2,			/* B, LDB: E */
		    one,			/* BETA:   1 */
		    Xx + ENTRY_SIZE*k1, d,      /* C, LDC: x1 */
                    Common->blas_ok) ;
	    }

	    /* solve L1'*x1 */
	    SUITESPARSE_BLAS_dtrsm ("L", "L", "C", "N",
		nscol,	nrhs,			/* M, N: x1 is nscol-by-nrhs */
		one,				/* ALPHA:  1 */
		Lx + ENTRY_SIZE*psx, nsrow,	/* A, LDA: L1 */
		Xx + ENTRY_SIZE*k1, d,          /* B, LDB: x1 */
                Common->blas_ok) ;

#else

	    /* x1 = x1 - L2'*E */
	    if (nsrow2 > 0)
	    {
		SUITESPARSE_BLAS_zgemm ("C", "N",
		    nscol, nrhs, nsrow2,	/* M, N, K */
		    minus_one,			/* ALPHA:  -1 */
		    Lx + ENTRY_SIZE*(psx + nscol),  /* A, LDA: L2 */
		    nsrow,
		    Ex, nsrow2,			/* B, LDB: E */
		    one,			/* BETA:   1 */
		    Xx + ENTRY_SIZE*k1, d,      /* C, LDC: x1 */
                    Common->blas_ok) ;
	    }

	    /* solve L1'*x1 */
	    SUITESPARSE_BLAS_ztrsm ("L", "L", "C", "N",
		nscol,	nrhs,			/* M, N: x1 is nscol-by-nrhs */
		one,				/* ALPHA:  1 */
		Lx + ENTRY_SIZE*psx, nsrow,	/* A, LDA: L1 */
		Xx + ENTRY_SIZE*k1, d,          /* B, LDB: x1 */
                Common->blas_ok) ;

#endif

	}
    }
}

#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
