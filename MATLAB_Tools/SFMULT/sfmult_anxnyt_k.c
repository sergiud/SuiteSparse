//==============================================================================
//=== sfmult_anxnyt_k ==========================================================
//==============================================================================

// SFMULT, Copyright (c) 2009, Timothy A Davis. All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

// y = (A*x)'	    where x has 2, 3, or 4 columns

// compare with sfmult_anxtyt_k

// sfmult_AN_XN_YT_2	y = (A*x)'  where x is n-by-2, and y is 2-by-m
// sfmult_AN_XN_YT_3	y = (A*x)'  where x is n-by-3, and y is 3-by-m (ldy = 4)
// sfmult_AN_XN_YT_4	y = (A*x)'  where x is n-by-4, and y is 4-by-m

#include "sfmult.h"

void sfmult_AN_XN_YT_2	// y = (A*x)'	x is n-by-2, and y is 2-by-m
(
    // --- outputs, not initialized on input
    double *Yx,		// 2-by-m
    double *Yz,		// 2-by-m if Y is complex (TO DO)

    // --- inputs, not modified
    const Int *Ap,	// size n+1 column pointers
    const Int *Ai,	// size nz = Ap[n] row indices
    const double *Ax,	// size nz values
    const double *Az,	// size nz imaginary values if A is complex (TO DO)
    Int m,		// A is m-by-n
    Int n,
    const double *Xx,	// n-by-2
    const double *Xz,	// n-by-2 if X complex (TO DO)
    int ac,		// true: use conj(A), otherwise use A (TO DO)
    int xc,		// true: use conj(X), otherwise use X (TO DO)
    int yc		// true: compute conj(Y), otherwise compute Y (TO DO)
)
{
    double x [2], a [2] ;
    Int p, pend, j, i0, i1 ;

    for (i0 = 0 ; i0 < m ; i0++)
    {
	Yx [2*i0  ] = 0 ;
	Yx [2*i0+1] = 0 ;
    }
    p = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	pend = Ap [j+1] ;
	x [0] = Xx [j  ] ;
	x [1] = Xx [j+n] ;
	if ((pend - p) % 2)
	{
	    i0 = Ai [p] ;
	    a [0] = Ax [p] ;
	    Yx [2*i0  ] += a [0] * x [0] ;
	    Yx [2*i0+1] += a [0] * x [1] ;
	    p++ ;
	}
	for ( ; p < pend ; p += 2)
	{
	    i0 = Ai [p  ] ;
	    i1 = Ai [p+1] ;
	    a [0] = Ax [p  ] ;
	    a [1] = Ax [p+1] ;
	    Yx [2*i0  ] += a [0] * x [0] ;
	    Yx [2*i0+1] += a [0] * x [1] ;
	    Yx [2*i1  ] += a [1] * x [0] ;
	    Yx [2*i1+1] += a [1] * x [1] ;
	}
    }
}


//==============================================================================
//=== sfmult_AN_XN_YT_3 ========================================================
//==============================================================================

void sfmult_AN_XN_YT_3	// y = (A*x)'	x is n-by-3, and y is 3-by-m (ldy = 4)
(
    // --- outputs, not initialized on input
    double *Yx,		// 3-by-m
    double *Yz,		// 3-by-m if Y is complex (TO DO)

    // --- inputs, not modified
    const Int *Ap,	// size n+1 column pointers
    const Int *Ai,	// size nz = Ap[n] row indices
    const double *Ax,	// size nz values
    const double *Az,	// size nz imaginary values if A is complex (TO DO)
    Int m,		// A is m-by-n
    Int n,
    const double *Xx,	// n-by-3
    const double *Xz,	// n-by-3 if X complex (TO DO)
    int ac,		// true: use conj(A), otherwise use A (TO DO)
    int xc,		// true: use conj(X), otherwise use X (TO DO)
    int yc		// true: compute conj(Y), otherwise compute Y (TO DO)
)
{
    double x [4], a [2] ;
    Int p, pend, j, i0, i1 ;

    for (i0 = 0 ; i0 < m ; i0++)
    {
	Yx [4*i0  ] = 0 ;
	Yx [4*i0+1] = 0 ;
	Yx [4*i0+2] = 0 ;
    }
    p = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	pend = Ap [j+1] ;
	x [0] = Xx [j    ] ;
	x [1] = Xx [j+  n] ;
	x [2] = Xx [j+2*n] ;
	if ((pend - p) % 2)
	{
	    i0 = Ai [p] ;
	    a [0] = Ax [p] ;
	    Yx [4*i0  ] += a [0] * x [0] ;
	    Yx [4*i0+1] += a [0] * x [1] ;
	    Yx [4*i0+2] += a [0] * x [2] ;
	    p++ ;
	}
	for ( ; p < pend ; p += 2)
	{
	    i0 = Ai [p  ] ;
	    i1 = Ai [p+1] ;
	    a [0] = Ax [p  ] ;
	    a [1] = Ax [p+1] ;
	    Yx [4*i0  ] += a [0] * x [0] ;
	    Yx [4*i0+1] += a [0] * x [1] ;
	    Yx [4*i0+2] += a [0] * x [2] ;
	    Yx [4*i1  ] += a [1] * x [0] ;
	    Yx [4*i1+1] += a [1] * x [1] ;
	    Yx [4*i1+2] += a [1] * x [2] ;
	}
    }
}


//==============================================================================
//=== sfmult_AN_XN_YT_4 ========================================================
//==============================================================================

void sfmult_AN_XN_YT_4	// y = (A*x)'	x is n-by-4, and y is 4-by-m
(
    // --- outputs, not initialized on input
    double *Yx,		// 4-by-m
    double *Yz,		// 4-by-m if Y is complex (TO DO)

    // --- inputs, not modified
    const Int *Ap,	// size n+1 column pointers
    const Int *Ai,	// size nz = Ap[n] row indices
    const double *Ax,	// size nz values
    const double *Az,	// size nz imaginary values if A is complex (TO DO)
    Int m,		// A is m-by-n
    Int n,
    const double *Xx,	// n-by-4
    const double *Xz,	// n-by-4 if X complex (TO DO)
    int ac,		// true: use conj(A), otherwise use A (TO DO)
    int xc,		// true: use conj(X), otherwise use X (TO DO)
    int yc		// true: compute conj(Y), otherwise compute Y (TO DO)
)
{
    double x [4], a [2] ;
    Int p, pend, j, i0, i1 ;

    for (i0 = 0 ; i0 < m ; i0++)
    {
	Yx [4*i0  ] = 0 ;
	Yx [4*i0+1] = 0 ;
	Yx [4*i0+2] = 0 ;
	Yx [4*i0+3] = 0 ;
    }
    p = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	pend = Ap [j+1] ;
	x [0] = Xx [j    ] ;
	x [1] = Xx [j+  n] ;
	x [2] = Xx [j+2*n] ;
	x [3] = Xx [j+3*n] ;
	if ((pend - p) % 2)
	{
	    i0 = Ai [p] ;
	    a [0] = Ax [p] ;
	    Yx [4*i0  ] += a [0] * x [0] ;
	    Yx [4*i0+1] += a [0] * x [1] ;
	    Yx [4*i0+2] += a [0] * x [2] ;
	    Yx [4*i0+3] += a [0] * x [3] ;
	    p++ ;
	}
	for ( ; p < pend ; p += 2)
	{
	    i0 = Ai [p  ] ;
	    i1 = Ai [p+1] ;
	    a [0] = Ax [p  ] ;
	    a [1] = Ax [p+1] ;
	    Yx [4*i0  ] += a [0] * x [0] ;
	    Yx [4*i0+1] += a [0] * x [1] ;
	    Yx [4*i0+2] += a [0] * x [2] ;
	    Yx [4*i0+3] += a [0] * x [3] ;
	    Yx [4*i1  ] += a [1] * x [0] ;
	    Yx [4*i1+1] += a [1] * x [1] ;
	    Yx [4*i1+2] += a [1] * x [2] ;
	    Yx [4*i1+3] += a [1] * x [3] ;
	}
    }
}
