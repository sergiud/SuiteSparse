/* ========================================================================== */
/* === PIRO_BAND/Source/piro_band_givens.c ================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for 
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Functions to find the Givens rotations to reduce a single entry to zero.
 * Both the real and complex versions uses LAPACK 3.0's algorithms to compute
 * the plane rotations reliably.
 * */

/* Global variables for the safe mininum (piro_band_safemin2 and 
 * piro_band_safemin2_f) and safe maximum (piro_band_safemx2 and
 * piro_band_safemx2_f) are expected to be set correctly before calling
 * PIRO_BAND(givens). */


#ifndef PIROBAND_COMPLEX /* [ */

/* Find the Givens rotations such that
 *  |  c    s |   | a  |   | d | 
 *  | -s    c | * | b  | = | 0 | 
 *  for a given a and b. 
 *  */
void PIRO_BAND(givens)
(
    Entry *a,               /* Entry array of size 1, i/p to rotate */
    Entry *b,               /* Entry array of size 1, i/p to rotate */
    Entry *c_op,            /* cosine of the givens rotation */
    Entry *s_op,            /* sine of the givens rotation */
    Entry *r_op,            /* rotated vector's first entry */
    PIRO_BAND_Common *Common /* common data structure */
)
{
    Entry f, g ;
    Entry c, s, r ;
    Entry f1, g1 ;
    Entry scale ;
    Int i, count ;
    Entry safemin2, safemx2 ;

    f = a[0] ;
    g = b[0] ;
    safemin2 = Common->safemin2 ;
    safemx2 = Common->safemx2 ;
    
    if (g == 0.0)
    {
        c = 1.0 ;
        s = 0.0 ;
        r = f ;
    }
    else if (f == 0.0)
    {
        c = 0.0 ;
        s = 1.0 ;
        r = g ;
    }
    else
    {
        f1 = f ;
        g1 = g ;
        scale = MAX (ABS(f1), ABS(g1)) ;
        count = 0 ;
        if (scale >= safemx2)
        {
            /* Avoid Overflow */
            while (scale >= safemx2)
            {
                count++ ;
                f1 = f1 * safemin2 ;
                g1 = g1 * safemin2 ;
                scale = MAX (ABS(f1), ABS(g1)) ;
            }
            r = sqrt (f1 * f1 + g1 * g1) ;
            c = f1 / r ;
            s = g1 / r ;
            for (i = 0 ; i < count ; i++)
            {
                r = r * safemx2 ;
            }
        }
        else if (scale <= safemin2)
        {
            /* Avoid Underflow */
            while (scale <= safemin2)
            {
                count++ ;
                f1 = f1 * safemx2 ;
                g1 = g1 * safemx2 ;
                scale = MAX (ABS(f1), ABS(g1)) ;
            }
            r = sqrt (f1 * f1 + g1 * g1) ;
            c = f1 / r ;
            s = g1 / r ;
            for (i = 0 ; i < count ; i++)
            {
                r = r * safemin2 ;
            }
        }
        else
        {
            /* Normal case */
            r = sqrt (f1 * f1 + g1 * g1) ;
            c = f1 / r ;
            s = g1 / r ;
        }

        /* Set the sign of c, s and r */
        if (ABS(f) > ABS(g) && c < 0)
        {
            c = -c ;
            s = -s ;
            r = -r ;
        }
    }

    /* Set the retrun values */
    c_op[0] = c ;
    s_op[0] = s ;
    r_op[0] = r ;
}


#else /* ] [ */

static Entry PIRO_BAND(abs1)
(
    Entry *ff,
    Int index
)
{
    Entry ff_real, ff_imag ;
    ff_real = ABS(ff[index]) ;
    ff_imag = ABS(ff[index+1]) ;
    return MAX(ff_real, ff_imag) ;
}

#define ABSSQ(f) (f[0] * f[0] + f[1] * f[1]) ;

/* Finding the hypotenuse accurately without the overflow.
 *
 * Algorithm 312 : Absolute value and Square root of a complex number, Paul
 * Friedland, Commun. of the ACM, Vol 10, No, 10, 1967.
 * */
Entry PIRO_BAND(hypot)
(
    Entry x,      
    Entry y
)
{
    Entry s, r ;

    x = ABS (x) ;
    y = ABS (y) ;
    if (x >= y)
    {
        if (x + y == x)
        {
            s = x ;
        }
        else
        {
            r = y / x ;
            s = x * sqrt (1.0 + r*r) ;
        }
    }
    else
    {
        ASSERT (y != 0.0) ;
        r = x / y ;
        s = y * sqrt (1.0 + r*r) ;
    } 
    return (s) ;
}

/* Find the Givens rotations such that
 *  |  c    s |   | a  |   | d | 
 *  | -s'   c | * | b  | = | 0 | 
 *  for a given a and b. 
 *
 *  c will always be real. But for a simpler indexing in storing and accessing
 *  the rotations we will use an array of size 2 with the imaginary part as
 *  zero.
 *  */
void PIRO_BAND(givens)
(
    Entry *farray,          /* Entry array of size 2, i/p to rotate */
    Entry *garray,          /* Entry array of size 2, i/p to rotate */
    Entry *c,               /* cosine of the givens rotation */
    Entry *s,               /* sine of the givens rotation */
    Entry *r,               /* rotated vector's first entry */
    PIRO_BAND_Common *Common /* common data structure */
)
{
    Entry f[2], g[2], ff[2], temp[2], temp2[2], r1[2] ;
    Entry d, di, dr, f2, f2s, g2, g2s ;
    Entry scale ;
    Entry af, ag ;
    Int k, count ;
    Entry safemin2, safemx2 ;

    safemin2 = Common->safemin2 ;
    safemx2 = Common->safemx2 ;
    ASSIGN_TO_SCALAR(f, farray, 0) ;
    ASSIGN_TO_SCALAR(g, garray, 0) ;

    count = 0 ;
    c[1] = 0.0 ; /* c is always real */

    if (g[0] == 0.0 && g[1] == 0.0)
    {
        /* Quick return */
        c[0] = 1.0 ;
        ASSIGN_ZERO_TO_MATRIX(s, 0) ;
        ASSIGN_TO_MATRIX(f, r, 0) ;
        return ;
    }

    af = PIRO_BAND(abs1)(f, 0) ;
    ag = PIRO_BAND(abs1)(g, 0) ;
    scale = MAX(af, ag) ;

    if (scale >= safemx2)
    {
        /* Avoid Overflow */
        while (scale >= safemx2)
        {
            count++ ;
            SCALE(f,  safemin2) ;
            SCALE(g,  safemin2) ;
            scale = scale * safemin2 ;
        }
    }
    else if (scale <= safemin2)
    {
        /* Avoid Underflow */
        while (scale <= safemin2)
        {
            count-- ;
            SCALE(f,  safemx2) ;
            SCALE(g,  safemx2) ;
            scale = scale * safemx2 ;
        }
    }

    f2 = ABSSQ(f) ;
    g2 = ABSSQ(g) ;

    if (f2 <= ((MAX(g2, 1)) * PIRO_BAND_Entry_MIN) )
    {
        /* f is very small */
        if (farray[0] == 0.0 && farray[1] == 0.0)
        {
            /* f is zero, quick return */
            /* c = 0, r = abs(garray) */
            c[0] = 0.0 ;
            r[0] = PIRO_BAND(hypot)(garray[0], garray[1]) ;
            r[1] = 0.0 ;
            d = PIRO_BAND(hypot)(g[0], g[1]) ;
            s[0] = g[0] / d ;
            s[1] = -g[1] / d ;
            return ;
        }
        /* c = abs(f)/abs(g) */
        f2s = PIRO_BAND(hypot)(f[0], f[1]) ;
        g2s = sqrt (g2) ;
        c[0] = f2s / g2s ;

        if (PIRO_BAND(abs1)(farray, 0) > 1)
        {
            d = PIRO_BAND(hypot)(farray[0], farray[1]) ;
            ff[0] = farray[0] /d ;
            ff[1] = farray[1] /d ;
        }
        else
        {
            dr = safemx2 * farray[0] ;
            di = safemx2 * farray[1] ;
            d = PIRO_BAND(hypot) (dr, di) ;
            ff[0] = dr / d ;
            ff[1] = di / d ;
        }
        /* s = ff * (g/g2s) */
        temp[0] = g[0]/g2s ;
        temp[1] = - g[1]/g2s ;
        MULT(s, ff, temp) ;

        /* r = c * d + s * g */
        MULT(r1, s, garray) ; 
        r[0] = r1[0] + c[0] * farray[0] ;
        r[1] = r1[1] + c[0] * farray[1] ;
    }
    else
    {
        /* Normal case */
        f2s = sqrt ( 1 + g2/f2) ;
        r1[0] = f2s * f[0] ;
        r1[1] = f2s * f[1] ; 

        c[0] = 1 / f2s ;

        d = f2 + g2 ;
        temp2[0] = r1[0]/d ;
        temp2[1] = r1[1]/d ;
        /* sn = sn * CONJ(g) ; */
        CONJ(temp, g) ;
        MULT(s, temp, temp2) ; 

        if ( count != 0)
        {
            if (count > 0)
            {
                for ( k = 0 ; k < count ; k++)
                {
                    SCALE(r1, safemx2) ;
                }
            }
            else
            {
                for ( k = 0 ; k < -count ; k++)
                {
                    SCALE(r1, safemin2) ;
                }
            }
        }
        ASSIGN_TO_MATRIX(r1, r, 0) ;
    }
    
}

#endif /* ] */
