/* ========================================================================== */
/* === PIRO_BAND/Demo/piro_band_testutils.c ================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * PIRO_BAND.  Version 1.0.
 * Copyright (C) 2009-2012, Sivasankaran Rajamanickam, Timothy A. Davis.
 * PIRO_BAND is licensed under Version 2.1 of the GNU Lesser
 * General Public License.  See lesser.txt for a text of the license.
 * PIRO_BAND is also available under other licenses; contact authors for 
 * details. http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Utility functions for testing the band reduction package. All the functions
 * assume C99 complex. 
 * */

#include "test_rand.h" 

/* Multiply the output from the band reduction, A1 = U * B * VT' */
void PIRO_BAND(svdmult)(
    COV_Entry *U, 
    COV_Entry *VT, 
    Entry *B1, 
    Entry *B2, 
    Int m, 
    Int n, 
    COV_Entry *A1, 
    Int sym
)
{
    COV_Entry *U1 ;
    int i, j, j2 , minmn ;
    int i1, i2, i3, i4 ;

    /* Allocate temproary space */
    U1 = (COV_Entry *) calloc(m*n, sizeof(COV_Entry)) ;
    if (!U1)
    {
        printf("Out of memory \n") ;
        return ; /* return failure : TBD */
    }

    minmn = m < n ? m : n ;

    /* multiply first column of U with B1[0] */
    for (i = 0 ; i < m ; i++)
    {
        if (!sym)
        {
            U1[i] = U[i] * ((COV_Entry) B1[0]) ;
        }
        else
        {
            U1[i] = U[i] * ((COV_Entry) B1[0]) + U[m+i] * ((COV_Entry) B2[1]) ;
        }
    }

    for (j = 1 ; j < minmn ; j++)
    {
        /* multiply columns j-1 and j of U with B2[k] and B1[k] and add
         *  to form column j */
        i1 = j ;
        i2 = j * m ;
        i3 = (j-1) * m ;
        i4 = (j+1) * m ;
        for (i = 0 ; i < m ; i++)
        {
            if (!sym || (sym && j == minmn-1) )
            {
                U1[i2+i] = U[i3+i]*((COV_Entry) B2[i1]) +
                                U[i2+i]*((COV_Entry) B1[j]) ;
            }
            else
            {
                U1[i2+i] = U[i3+i]*((COV_Entry) B2[i1]) +
                            U[i2+i]*((COV_Entry) B1[j]) +
                            U[i4+i]*((COV_Entry) B2[i1+1]) ;
            }
        }
    }

    for (j = 0 ; j < n ; j++)
    {
        for (i = 0 ; i < m ; i++)
        {
            A1[j*m+i] = U1[i] * VT[j*n] ;
        }
        for (j2 = 1 ; j2 < n ; j2++)
        {
            for (i = 0 ; i < m ; i++)
            {
                A1[j*m+i] += U1[j2*m+i] * VT[j*n+j2] ;
            }
        }
    }

    /* Free temporary space */
    free(U1) ;
}

/* Find the column norm of a dense matrix */
Entry PIRO_BAND(find_norm)(COV_Entry *A, Int m, Int n)
{
    Int i, j ;
    Entry norm, current ;

    norm = 0 ;
    for (j = 0 ; j < n ; j++)
    {
        current = 0.0 ;
        for (i = 0 ; i < m ; i++)
        {
            current += COV_ABS(A[i+j*m]) ;
        }
        if (current > norm) /* NAN ?? */
        {
            norm = current ;
        }
    }
    return norm ;
}

/* Find the column norm of a band matrix */
Entry PIRO_BAND(find_band_norm)(COV_Entry *A, Int ldab, Int m, Int n, Int bl,
                                Int bu)
{
    Int i, j ;
    Int start, end ;
    Entry norm, current ;

    norm = 0 ;
    for (j = 0 ; j < n ; j++)
    {
        current = 0.0 ;
        start = MAX(j - bu, 0) ;
        end = MIN(j + bl, m-1) ;
        for (i = start ; i <= end ; i++)
        {
            current += COV_ABS(A[INDEX(i, j)]) ;
        }
        if (current > norm) /* NAN ?? */
        {
            norm = current ;
        }
    }
    return norm ;
}

/* Read a band matrix stored in a file in packed format */
Int PIRO_BAND(get_matrix)(char *pfile, Int ldab, Int n, COV_Entry *A)
{
    float dtemp = 0 ;
#ifdef PIROBAND_COMPLEX
    float dtemp1 = 0 ;
#endif
    FILE *fp1 ;
    Int i1, j, ok = 1 ;

    fp1 = fopen(pfile, "r") ;
    if (!fp1)
    {
        printf("File %s not found\n", pfile) ;
        return 0 ;
    }
    for (j = 0 ; j < n && ok ; j++)
    {
        for (i1 = 0 ; i1 < ldab && ok ; i1++)
        {
#ifdef PIROBAND_COMPLEX
            ok = fscanf(fp1, "%f %f", &dtemp, &dtemp1) == 2 ;
            A[i1+(j*ldab)] = (COV_Entry) (dtemp + dtemp1*I) ;
#else
            ok = fscanf(fp1, "%f", &dtemp) == 1 ;
            A[i1+(j*ldab)] = (COV_Entry) dtemp ;
#endif
        }
    }
    if (!ok)
    {
        printf ("invalid file %s\n", pfile) ;
    }
    fclose(fp1) ;
    return ok ;
}

/* Generate a random band matrix in packed format */
Int PIRO_BAND(get_random_matrix)(Int ldab, Int n, COV_Entry *A, Int sym, Int bu)
{
    Int i1, j ;

    for (j = 0 ; j < n ; j++)
    {
        for (i1 = 0 ; i1 < ldab ; i1++)
        {
#ifdef PIROBAND_COMPLEX
            if (sym && i1 == bu)
            {
                A[i1+(j*ldab)] = (COV_Entry) (xrand(1.0) + 0.0 * I) ;
            }
            else
            {
                A[i1+(j*ldab)] = (COV_Entry) (xrand(1.0) + xrand(1.0) * I) ;
            }
#else
            A[i1+(j*ldab)] = (COV_Entry) xrand(1.0) ;
#endif
        }
    }
    return 1 ;
}

/* Allocate space for U, V, and C to accumulate the Givens rotations and
 * initialize c to identity
 * */
Int PIRO_BAND(init_input)(
    Int m,
    Int n,
    COV_Entry **temp1,
    COV_Entry **U,
    COV_Entry **V,
    COV_Entry **C
)
{
    COV_Entry *temp ;
    *temp1 = (COV_Entry *) malloc (((2 * m * m) + (n * n)) * sizeof(COV_Entry));
    temp = *temp1 ;
    if (!temp)
    {
        printf("Out of memory") ;
        return 0 ;
    }
    *U = temp ;
    temp += (m*m) ;
    *C = temp ;
    temp += (m*m) ;
    *V = temp ;
    temp += (n*n) ;
    PIRO_BAND(identity_matrix) (m, *C) ;
    return 1 ;
}

void PIRO_BAND(identity_matrix)
(
    Int m,
    COV_Entry *C
)
{
    /* Initialize C to identity */
    Int i, j ;
    for (j = 0 ; j < m ; j++)
    {
        C [j*m+j] = 1.0 + 0.0*I;
        for (i = 0 ; i < m ; i++)
        {
            if (i != j) C [j*m+i] = 0.0 + 0.0*I ;
        }
    }
}

/* Check the result after band reduction.
 * If C is not NULL it should equal U' or U.' on input.  It is overwritten with
 * C-U' or C-U.'.  If lapack is false then V is overwritten with V' */
Int PIRO_BAND(chk_output)
(
    Int m,          /* A is m-by-n */
    Int n,
    Int bu,         /* upper bandwidth of A */
    Int bl,         /* lower bandwidth of A */
    Int ldab,       /* leading dimension of A */
    Int sym,        /* 0: A is unsymmetric, 1: A is symmetric */
    Entry *D,       /* diagonal of bidiagonal form, B */
    Entry *E,       /* upper diagonal of bidiagonal form */
    COV_Entry *U,   /* A = U*B*V' should hold */
    COV_Entry *V,
    COV_Entry *C,   /* piro_band_reduce reductions also applied to C */
    COV_Entry *A,   /* A, as provided on input to piro_band_reduce */
    Int lapack      /* 0: C=U.', 1: C=U', do not check C if C is NULL */
)
{
    Int i, j, start, end ;
    COV_Entry *A1 ;
    Entry norm ;
    Entry anorm ;

    if (C != NULL)
    {
        /* We use identity for C on input to piro_band_reduce,
         * so on output from piro_band_reduce, C should be the same as U after
         * reduction and the transpose of U in the lapack interfaces.
         * Compute C = C-U' or C-U.' */
        for (j = 0 ; j < m ; j++)
        {
            if (lapack)
            {
                C[j*m+j] -= COV_CONJ(U[j*m+j]) ;
            }
            else
            {
                C[j*m+j] -= U[j*m+j] ;
            }
            for ( i = 0 ; i < m ; i++)
            {
                if (i != j)
                {
                    if (lapack)
                    {
                        C[i*m+j] -= COV_CONJ(U[j*m+i]) ;
                    }
                    else
                    {
                        C[i*m+j] -= U[i*m+j] ;
                    }
                }
            }
        }

        /* Find the norm of the difference */
        printf("C - U' ----") ;
        norm = PIRO_BAND(find_norm)(C, m, m) ;
        printf("Norm is %0.6g\n", norm) ;
        if ((sizeof (Entry) > 4) ? (norm > 1e-12) : (norm > 1e-5))
        {
            printf("Norm is too high!\n") ;
            return 0 ;
        }
    }

    /* Allocate temp space for the result */
    A1 = (COV_Entry *) calloc(m * n, sizeof(COV_Entry)) ;
    if (!A1)
    {
        printf("Out of memory") ;
        return 0 ;
    }
    if (!lapack)
    {
        /* Transpose V for the reduce */
        PIRO_BAND(inplace_conjugate_transpose)(n, (Entry *) V, n) ;
    }

    /* Multiply the r.h.s A1 = U * B * V' */
    PIRO_BAND(svdmult)(U, V, D, E, m, n, A1, sym) ;

    /* Find the difference between the copy of the band matrix and A1 */
    for (j = 0 ; j < n ; j++)
    {
        start = MAX(j - bu, 0) ;
        end = MIN(j + bl, m-1) ;
        for (i = start ; i <= end ; i++)
        {
            A1[i+(j*m)] = A[INDEX(i, j)] - A1[i+(j*m)] ;
            if (sym && i != j)
            {
                A1[j+(i*m)] = A[INDEX(i, j)] - conj(A1[j+(i*m)]) ;
            }
        }
    }

    /* Find the norm of the difference */
    printf("A - U * B * V' ----") ;
    norm = PIRO_BAND(find_norm)(A1, m, n) ;
    anorm = PIRO_BAND(find_band_norm)(A, ldab, m, n, bl, bu) ;
    norm = norm / anorm ;
    printf("Norm is %0.6g\n",norm) ;
    free(A1) ;
    if ((sizeof (Entry) > 4) ? (norm > 1e-12) : (norm > 1e-5))
    {
        printf("Norm is too high!\n") ;
        return 0 ;
    }
    return 1 ;
}

/* Print the output from band reductio */
void PIRO_BAND(print_output)
(
    Int m,
    Int n,
    Entry *D,
    Entry *E,
    COV_Entry *U,
    COV_Entry *V,
    COV_Entry *C
)
{
    Int i ;

    for ( i = 0 ; i < MIN(m, n) ; i++ )
    {
        printf("D[%g] = %0.8f\n", (double) i, D[i]) ;
    }
    for ( i = 0 ; i < MIN(m, n) ; i++ )
    {
        printf("E[%g] = %0.8f\n", (double) i, E[i]) ;
    }
#ifdef PIROBAND_COMPLEX
    for ( i = 0 ; i < m*m ; i++ )
    {
        printf("U[%g] = %0.8f %0.8f\n", (double) i, creal(U[i]), cimag(U[i])) ;
    }
    for ( i = 0 ; i < m*m ; i++ )
    {
        printf("C[%g] = %0.8f %0.8f\n", (double) i, creal(C[i]), cimag(C[i])) ;
    }
    for ( i = 0 ; i < n*n ; i++ )
    {
        printf("V[%g] = %0.8f %0.8f\n", (double) i, creal(V[i]), cimag(V[i])) ;
    }
#else
    for ( i = 0 ; i < m*m ; i++ )
    {
        printf("U[%g] = %0.8f\n", (double) i, U[i]) ;
    }
    if (C != NULL)
    {
        for ( i = 0 ; i < m*m ; i++ )
        {
            printf("C[%g] = %0.8f\n", (double) i, C[i]) ;
        }
    }
    for ( i = 0 ; i < n*n ; i++ )
    {
        printf("V[%g] = %0.8f\n", (double) i, V[i]) ;
    }
#endif
}
