/* ========================== genband_givens.c ============================= */
/* Functions to find the Givens rotations to reduce a single entry to zero.
 * Both the real and complex versions are based on LAPACK's dlartg.f and
 * clartg.f algorithms.

 TODO use piro_band instead?
 * */

#include "genband_util.h"
#include "genband_internal.h"
#include "blksky.h"

#ifndef PIROBAND_COMPLEX /* [ */

void apply_givens
(
    double *da,
    double *db,
    double c,
    double s
)
{
    double t_da, t_db, temp,  temp2 ;

    t_da = da[0] ;
    t_db = db[0] ;

    if (c == 0.0 && s == 1.0)
    {
        /* swap */
        da[0] = t_db ;
        db[0] = t_da ;
        return ;
    }

    if (t_da == 0.0)
    {
        if (t_db == 0.0)
        {
            return ;
        }
        else
        {
            t_da = t_db * s ; 
            t_db = t_db * c ;
        }
    }
    else if (t_db == 0.0)
    {
        t_db = -t_da * s ;
        t_da = t_da * c ;
    }
    else
    {
        temp = t_da ;
        temp2 = t_db ;
        t_da = temp * c ;
        t_da += t_db * s ; 
        t_db = temp2 * c ;
        t_db -= temp * s ;
    }

    da[0] = t_da ;
    db[0] = t_db ;
    return ;
}

void givens
(
    double *a,
    Int a_index,
    double *b,    
    Int b_index,
    double *c_ip,
    Int c_index,
    double *s_ip,
    Int s_index,
    double *r_op,
    Int r_index,
    sky_common *scom
)
{
    /*double safemin, safemin2, safemx2, eps ;*/
    double safemin2, safemx2 ;
    double f, g ;
    double c, s, r ;
    double f1, g1 ;
    double scale ;
    /* Int esfmn2 ;*/
    Int i, count ;

    f = a [a_index] ;
    g = b [b_index] ;
    
#if 0
    /* TBD : Should adjust for float */
    safemin = GENBAND_DBL_MIN ;
    eps = GENBAND_EPS/2 ;
    esfmn2  = log(safemin/eps) / log(2) / 2 ;
    safemin2 = pow(2, esfmn2) ; /* TBD : powf and -ansi won't work together */
    safemx2 = 1/safemin2 ;
#endif

    safemin2 = scom->safemin2 ;
    safemx2 = scom->safemx2 ;

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
        if (scale >= safemx2)
        {
            count = 0 ;
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
            count = 0 ;
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
            r = sqrt (f1 * f1 + g1 * g1) ;
            c = f1 / r ;
            s = g1 / r ;
        }
        if (ABS(f) > ABS(g) && c < 0)
        {
            c = -c ;
            s = -s ;
            r = -r ;
        }
    }

    c_ip[c_index] = c ;
    s_ip[s_index] = s ;
    r_op[r_index] = r ;
}


#else /* ] [ */

static Entry GENBAND(abs1)
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

static Entry GENBAND(abssq)
(
    Entry *ff
)
{
    return ( ff[0] * ff[0] + ff[1] * ff[1] ) ;
}

Entry GENBAND(slapy2)
(
    Entry x,
    Entry y
)
{
    Entry xabs, yabs, w, z ;

    xabs = ABS(x) ;
    yabs = ABS(y) ;

    w = MAX(xabs, yabs) ;
    z = MIN(xabs, yabs) ;
    if ( z + w == w)
    {
        return w ;
    }
    else
    {
        return w * sqrt ( 1 + (z/w)*(z/w)) ; /* TBD : w != 0 */
    }
}

void GENBAND(givens)
(
    Entry *farray,
    Int f_index,
    Entry *garray,            /* usually same A, but leave it here */
    Int g_index,
    Entry *c,
    Int c_index,
    Entry *s,
    Int s_index,
    Entry *r,
    Int r_index
)
{
    Entry f[2], g[2], ff[2], temp[2], temp2[2] ;
    Entry safemin, safemin2, safemx2, eps ;
    Entry d, di, dr, f2, f2s, g2, g2s ;
    Entry scale ;
    Entry af, ag ;
    Entry *tmpptr ;
    Int esfmn2 ;
    Int k, count ;

    safemin = GENBAND_DBL_MIN ;
    eps = GENBAND_EPS/2 ;
    esfmn2  = log(safemin/eps) / log(2) / 2 ;
    safemin2 = pow(2, esfmn2) ; /* TBD : powf and -ansi won't work together */
    safemx2 = 1/safemin2 ;


    f[0] = farray [f_index] ;
    f[1] = farray [f_index+1] ;

    g[0] = garray [g_index] ;
    g[1] = garray [g_index+1] ;

    af = GENBAND(abs1)(f, 0) ;
    ag = GENBAND(abs1)(g, 0) ;
    scale = MAX(af, ag) ;

    count = 0 ;

    if (scale >= safemx2)
    {
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
        if (g[0] == 0.0 && g[1] == 0.0)
        {
            c[c_index] = 1 ;
            c[c_index+1] = 0 ;
            s[s_index] = 0 ;
            s[s_index+1] = 0 ;
            return ;
        }
        while (scale <= safemin2)
        {
            count-- ;
            SCALE(f,  safemx2) ;
            SCALE(g,  safemx2) ;
            scale = scale * safemx2 ;
        }
    }

    f2 = GENBAND(abssq) (f) ;
    g2 = GENBAND(abssq) (g) ;

    if (f2 < ((MAX(g2, 1)) * safemin) )
    {
        /* f is very small, rare case. */
        if (farray[f_index] == 0.0 && farray[f_index+1] == 0.0)
        {
            c[c_index] = 0.0 ;
            c[c_index+1] = 0.0 ;
            r[0] = GENBAND(slapy2)(garray[g_index], garray[g_index+1]) ;
            r[1] = 0.0 ;
            d = GENBAND(slapy2)(g[0], g[1]) ;
            s[s_index] = g[0] / d ;
            s[s_index+1] = -g[1] / d ;
            return ;
        }
        f2s = GENBAND(slapy2)(f[0], f[1]) ;
        g2s = sqrt (g2) ;
        c[c_index] = f2s / g2s ;
        c[c_index+1] = 0.0 ;

        if (GENBAND(abs1)(farray, f_index) > 1)
        {
            d = GENBAND(slapy2)(farray[f_index], farray[f_index+1]) ;
            ff[0] = farray[f_index] /d ;
            ff[1] = farray[f_index+1] /d ;
        }
        else
        {
            dr = safemx2 * farray[f_index] ;
            di = safemx2 * farray[f_index+1] ;
            d = GENBAND(slapy2) (dr, di) ;
            ff[0] = dr / d ;
            ff[1] = di / d ;
        }
        temp[0] = g[0]/g2s ;
        temp[1] = - g[1]/g2s ;
        tmpptr = s+s_index ;
        MULT(tmpptr, ff, temp) ;
        /* Fix the inices when moved to scalars ... */
        temp[0] = garray[g_index] ;
        temp[1] = garray[g_index+1] ;
        MULT(r, s, temp) ; 
        tmpptr = c+(c_index) ;
        temp[0] = farray[f_index] ;
        temp[1] = farray[f_index+1] ;
        MULT_ADD(r, tmpptr, temp) ; /* Assumes r has always 0, 1 index */
    }
    else
    {
        /* Normal case */
        f2s = sqrt ( 1 + g2/f2) ;
        r[0] = f2s * f[0] ;
        r[1] = f2s * f[1] ; 
        c[c_index] = 1 / f2s ;
        c[c_index+1] = 0.0 ;
        d = f2 + g2 ;
        s[s_index] = r[0]/d ;
        s[s_index+1] = r[1]/d ;
        /* sn = sn * CONJ(g) ; */
        CONJ(temp, g) ;
        temp2[0] = s[s_index] ;
        temp2[1] = s[s_index+1] ;
        tmpptr = s+s_index ;
        MULT(tmpptr, temp, temp2) ; 

        if ( count != 0)
        {
            if (count > 0)
            {
                for ( k = 0 ; k < count ; k++)
                {
                    SCALE(r, safemx2) ;
                }
            }
            else
            {
                for ( k = 0 ; k < -count ; k++)
                {
                    SCALE(r, safemin2) ;
                }
            }
        }
    }
    
}

#endif /* ] */
