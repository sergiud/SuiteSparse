//==============================================================================
//=== sfmult_AN_XT_YT ==========================================================
//==============================================================================

// SFMULT, Copyright (c) 2009, Timothy A Davis. All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

// y = (A*x')'	    A is m-by-n, x is k-by-n, y is k-by-m

// compare with sfmult_AN_XT_YN for kernel usage
// compare with sfmult_AN_XN_YT for outer loop structure but different kernels

#include "sfmult.h"

mxArray *sfmult_AN_XT_YT    // y = (A*x')'
(
    const mxArray *A,
    const mxArray *X,
    int ac,		// if true: conj(A)   if false: A. ignored if A real
    int xc,		// if true: conj(x)   if false: x. ignored if x real
    int yc		// if true: conj(y)   if false: y. ignored if y real
)
{
    mxArray *Y ;
    double *Ax, *Az, *Xx, *Xz, *Yx, *Yz, *Wx, *Wz ;
    Int *Ap, *Ai ;
    Int m, n, k, k1, i ;
    int Acomplex, Xcomplex, Ycomplex ;

    //--------------------------------------------------------------------------
    // get inputs
    //--------------------------------------------------------------------------

    m = mxGetM (A) ;
    n = mxGetN (A) ;
    k = mxGetM (X) ;
    if (n != mxGetN (X)) sfmult_invalid ( ) ;
    Acomplex = mxIsComplex (A) ;
    Xcomplex = mxIsComplex (X) ;
    Ap = mxGetJc (A) ;
    Ai = mxGetIr (A) ;
    Ax = mxGetPr (A) ;
    Az = mxGetPi (A) ;
    Xx = mxGetPr (X) ;
    Xz = mxGetPi (X) ;

    //--------------------------------------------------------------------------
    // allocate result
    //--------------------------------------------------------------------------

    Ycomplex = Acomplex || Xcomplex ;
    Y = sfmult_yalloc (k, m, Ycomplex) ;
    Yx = mxGetPr (Y) ;
    Yz = mxGetPi (Y) ;

    //--------------------------------------------------------------------------
    // special cases
    //--------------------------------------------------------------------------

    if (k == 0 || m == 0 || n == 0 || Ap [n] == 0)
    {
	// Y = 0
	return (sfmult_yzero (Y)) ;
    }
    if (k == 1)
    {
	// Y = A*X
	sfmult_AN_x_1 (Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc) ;
	return (Y) ;
    }
    if (k == 2)
    {
	// Y = (A * X')'
	sfmult_AN_XT_YT_2 (Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
	return (Y) ;
    }
    if (k == 4)
    {
	// Y = (A * X')'
	sfmult_AN_XT_YT_4 (Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
	return (Y) ;
    }

    if (k > 12 && Ap [n] < 8 * MAX (m,n))	// (TO DO) check overflow
    {
	// Y = (A * X')' when A is moderately sparse and X is large
	mxArray *C ;
	double *Cx, *Cz ;
	Int *Cp, *Ci ;
	// C = A' ;
	C = ssmult_transpose (A, 0) ;
	// (TO DO) if C is NULL, skip this and try using (A*X')' below
	// Y = (C' * X')' when A is moderately sparse and X is large
	Cp = mxGetJc (C) ;
	Ci = mxGetIr (C) ;
	Cx = mxGetPr (C) ;
	Cz = mxGetPi (C) ;
	sfmult_xA (Yx, Yz, Cp, Ci, Cx, Cz, n, m, Xx, Xz, ac, xc, yc, k) ;
	mxDestroyArray (C) ;
	return (Y) ;
    }

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    sfmult_walloc (4, m, &Wx, &Wz) ;
    // (TO DO) if walloc fails, use a workspace-free technique.
    // This may require new sparse-times-vector kernels (both x and y with
    // non-unit strides)

    //--------------------------------------------------------------------------
    // Y = (A*X')', in blocks of up to 4 columns of X, using sfmult_anxtyt
    //--------------------------------------------------------------------------

    k1 = k % 4 ;
    if (k1 == 1)
    {
	// W = A * X(1,:)'
	sfmult_AN_xk_1 (Wx, Wz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k) ;
	// Y (1,:) = W'
	for (i = 0 ; i < m ; i++)
	{
	    Yx [k*i] = Wx [i] ;
	}
	Yx += 1 ;
	Yz += 1 ;
	Xx += 1 ;
	Xz += 1 ;
    }
    else if (k1 == 2)
    {
	// W = (A * X(1:2,:)')'
	sfmult_AN_XT_YT_2 (Wx, Wz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
	// Y (1:2,:) = W
	for (i = 0 ; i < m ; i++)
	{
	    Yx [k*i  ] = Wx [2*i  ] ;
	    Yx [k*i+1] = Wx [2*i+1] ;
	}
	Yx += 2 ;
	Yz += 2 ;
	Xx += 2 ;
	Xz += 2 ;
    }
    else if (k1 == 3)
    {
	// W = (A * X(1:3,:)')'
	sfmult_AN_XT_YT_3 (Wx, Wz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
	// Y (1:3,:) = W
	for (i = 0 ; i < m ; i++)
	{
	    Yx [k*i  ] = Wx [4*i  ] ;
	    Yx [k*i+1] = Wx [4*i+1] ;
	    Yx [k*i+2] = Wx [4*i+2] ;
	}
	Yx += 3 ;
	Yz += 3 ;
	Xx += 3 ;
	Xz += 3 ;
    }
    for ( ; k1 < k ; k1 += 4)
    {
	// W = (A * X(1:4,:)')'
	sfmult_AN_XT_YT_4 (Wx, Wz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
	// Y (k1+(1:4),:) = W
	for (i = 0 ; i < m ; i++)
	{
	    Yx [k*i  ] = Wx [4*i  ] ;
	    Yx [k*i+1] = Wx [4*i+1] ;
	    Yx [k*i+2] = Wx [4*i+2] ;
	    Yx [k*i+3] = Wx [4*i+3] ;
	}
	Yx += 4 ;
	Yz += 4 ;
	Xx += 4 ;
	Xz += 4 ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    mxFree (Wx) ;
    return (Y) ;
}
