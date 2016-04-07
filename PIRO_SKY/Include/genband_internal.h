/* ================== genband_internal.h ==================================== */
/* 
 * Internal macros for general band reduction. All definitions in this file
 * depends on macros PIROBAND_LONG, PIROBAND_FLOAT, PIROBAND_COMPLEX defined
 * for int/long, double/float and real/complex support respectively.  All
 * function names have a suffix in the form _xyz to support the different
 * versions based on the defines. For eg, _dri will accept double, real and
 * integer arguments.  The real version of the all the macros will use an array
 * of size one and the complex version will use an array of size 2.  */

#undef ZERO
#undef ONE
#undef Int
#undef ID
#undef Double
#undef DID
#undef CONJ
#undef ABS
#undef GENBAND
#undef GENBAND_LONG_NAME
#undef GENBAND_LAPACK_NAME
#undef MULT
#undef MULT_ADD
#undef MULT_SUB
#undef GIVENS
#undef APPLY_GIVENS
#undef APPLY_GIVENS_TO_COL
#undef PRINT_NUMVALUES
#undef CSIZE
#undef ASSIGN_TO_SCALAR
#undef ASSIGN_TO_ZERO
#undef NEG
#undef ASSIGN_TO_MATRIX
#undef ASSIGN_ZERO_TO_MATRIX
#undef ASSIGN_MATRIX_TO_MATRIX
#undef ASSIGN_CS_SCALARS
#undef SAVE_CS_SCALARS
#undef SCALE
#undef INDEX
#undef AINDEX
#undef GENBAND_DBL_MIN
#undef GENBAND_EPS
#undef GENBAND_SAFEMIN2
#undef GENBAND_SAFEMX2  

#undef PIROBAND_LONG
#undef PIROBAND_FLOAT
#undef PIROBAND_COMPLEX


/* ========================================================================= */
/* ======= Definitions for int/long usage  ================================= */
/* ========================================================================= */

#include "SuiteSparse_config.h"

/* always use SuiteSparse_long for PIRO_SKY */
#define PIROBAND_LONG
#define Int SuiteSparse_long
#define ID "%ld"
#define GENBAND_LONG_NAME(name) genband_ ## name ## _l
#define GENBAND_LAPACK_NAME(a) genband_ ## a ## l

/* ========================================================================= */
/* ======= Definitions for double/float usage  ============================= */
/* ========================================================================= */

/* always use double version (not float) */

/*
#ifdef PIROBAND_FLOAT

#define Double float
#define DID "%0.4f"
#define GENBAND_DBL_MIN FLT_MIN
#define GENBAND_EPS FLT_EPSILON
#define GENBAND_SAFEMIN2 genband_safemin2_f
#define GENBAND_SAFEMX2 genband_safemx2_f

#else
*/

#define Double double
#define DID "%0.4f"
#define GENBAND_DBL_MIN DBL_MIN
#define GENBAND_EPS DBL_EPSILON
#define GENBAND_SAFEMIN2 genband_safemin2
#define GENBAND_SAFEMX2 genband_safemx2

/*
#endif
*/

#define ZERO 0.0
#define ONE 1.0

/* ========================================================================= */
/* ==================== Definitions for GENBAND(name) ====================== */
/* ========================================================================= */

/* The macro GENBAND(name)  gets defined to genband_name_xyz where 
 * x - s/d for single or double,  y - l/i for long or int, z - c/r for complex 
 * or real.
 * */

#ifdef PIROBAND_LONG

#ifdef PIROBAND_FLOAT

#ifdef PIROBAND_COMPLEX

#define GENBAND(name) genband_ ## name ## _scl

#else /* PIROBAND_COMPLEX */

#define GENBAND(name) genband_ ## name ## _srl

#endif /* PIROBAND_COMPLEX */

#else /* PIROBAND_FLOAT */

#ifdef PIROBAND_COMPLEX

#define GENBAND(name) genband_ ## name ## _dcl

#else /* PIROBAND_COMPLEX */

#define GENBAND(name) genband_ ## name ## _drl

#endif /* PIROBAND_COMPLEX */

#endif /* PIROBAND_FLOAT */

#else /* PIROBAND_LONG */

#ifdef PIROBAND_FLOAT

#ifdef PIROBAND_COMPLEX

#define GENBAND(name) genband_ ## name ## _sci

#else /* PIROBAND_COMPLEX */

#define GENBAND(name) genband_ ## name ## _sri

#endif /* PIROBAND_COMPLEX */

#else /* PIROBAND_FLOAT */

#ifdef PIROBAND_COMPLEX

#define GENBAND(name) genband_ ## name ## _dci

#else /* PIROBAND_COMPLEX */

#define GENBAND(name) genband_ ## name ## _dri

#endif /* PIROBAND_COMPLEX */

#endif /* PIROBAND_FLOAT */

#endif /* PIROBAND_LONG */


/* Entry can be one of the following datatype : double, float based on whether
 * PIROBAND_FLOAT and PIROBAND_COMPLEX are defined or not.
 * */
#define Entry  Double

int GENBAND(reduce)
(
    Int blks[],
    /*blk[0], blk[1],	       #columns, #rows in block for upper band */
    /*blk[2], blk[3],	       #columns, #rows in block for lower band */
    Int m,		    /* #rows in the original matrix            */
    Int n,		    /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,		    /* lower bandwidth                         */
    Int bu, 		    /* upper bandwidth                         */
    Entry *A,	            /* Band Matrix	                       */
    Int ldab,		    /* leading dimension of A                 */
    Entry *B1,              /* o/p diagonal 0                          */
    Entry *B2,              /* o/p diagonal 1                          */
    Entry *U,               /* o/p accumalated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumalated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Entry *dws,             /* workspace of size 2*MAX(nc*nr, ncl*nrl) */
    Int sym                 /* Is the matrix symmetric ?               */
) ;

/* ========================================================================= */
/* ======= Definitions to find Collumn and Row Givens' rotations ==========  */
/* ========================================================================= */

#ifndef PIROBAND_COMPLEX /* [ */
/*#define GIVENS(da, db, d, c, s) { GENBAND(givens) (da, db, c, s, d) ; \
                    PRINT(("gc da=%0.4f, db=%0.4f, c=%0.4f, s=%0.4f d=%0.4f\n",\
				da[0], db[0], c[0], s[0], d[0])) ;\
                        }*/
#define GIVENS(da, db, d, c, s) { \
                givens(&da, 0, &db, 0, &c, 0, &s, 0, &d, 0, scom ) ; \
                PRINT(("gc da=%0.4f, db=%0.4f,  c=%0.4f, s=%0.4f d=%0.4f\n",\
				da, db, c, s, d)) ;\
                        }


#else /* !PIROBAND_COMPLEX ][ */

#define GIVENS(da, db, d, c, s) {  GENBAND(givens) (da, db, c, s, d) ; \
PRINT(("gc da=%0.4f %0.4fi, db=%0.4f %0.4fi, c=%0.4f%0.4fi, s=%0.4f %0.4fi\n", \
                        da[0], da[1], db[0], db[1], c[0], c[1], s[0], s[1])) ;\
			}

#endif /* !PIROBAND_COMPLEX ] */

/* ========================================================================= */
/* ======= Definitions to apply Column and Row Givens' rotations ============ */
/* ========================================================================= */

/* Can be easily combined as CONJ handles REAL cases too. But leaving like this
 * for performance reasons as it will be one assignmnet less. c is a real. So a
 * complex multiply won't be necessary. Using temp2 in real case for performance
 * reasons.
 * */

#ifndef PIROBAND_COMPLEX /* [ */ 
/*#define APPLY_GIVENS(da, db, c, s, temp, conjs) { \
                        temp = da ; \
                        temp2 = db ; \
                        da = temp * c ; \
                        da += db * s ; \
                        db = temp2 * c ; \
                        db -= temp * s ; \
                        }*/

#define APPLY_GIVENS(da, db, c, s, temp, conjs) { \
                        apply_givens(&da, &db, c, s) ;\
                        }

/*#define APPLY_GIVENS(da, db, c, s, temp, conjs) { \
                        PRINT_NUMVALUES("b4 ac", da, db, c, s) ; \
			ASSIGN_TO_SCALAR(temp, da, 0) ; \
			ASSIGN_TO_SCALAR(temp2, db, 0) ; \
			MULT(da, temp, c) ; \
			MULT_ADD(da, db, s) ; \
			MULT(db, temp2, c) ; \
			MULT_SUB(db, temp, s) ; \
                        PRINT_NUMVALUES("ac", da, db, c, s) ; \
			}*/


#else /* !PIROBAND_COMPLEX ][ */

/*#define APPLY_GIVENS(da, db, c, s, temp, conjs) { \
                        PRINT_NUMVALUES("b4 ac", da, db, c, s) ; \
                        CONJ(conjs, s) ; \
			ASSIGN_TO_SCALAR(temp, da, 0) ; \
			SCALE(da, c[0]) ; \
			MULT_ADD(da, db, s) ; \
			SCALE(db, c[0]) ; \
			MULT_SUB(db, temp, conjs) ; \
                        PRINT_NUMVALUES("ac", da, db, c, s) ; \
			}*/

#endif /* !PIROBAND_COMPLEX ] */

/* Definition to Apply the Givens rotation to an entire column */
#define APPLY_GIVENS_TO_COL(U, ldu, col1, col2, si, ei, c, s, i, da, db, dtemp, n) { \
    ASSERT(si >=0 && si < n, "") ; \
    ASSERT(ei >=0 && ei < n, "") ; \
		/*mexPrintf("c=%g, s=%g\n ", c, s) ; \*/\
    for (i = si ; i <= ei ; i++) \
    {\
		da = U[AINDEX(i, col1, ldu)] ;\
		db = U[AINDEX(i, col2, ldu)] ;\
		/*mexPrintf("da=%g, db=%g ", da, db) ; \*/\
        APPLY_GIVENS(da, db, c, s, dtemp, conjs) ; \
		/*mexPrintf("da=%g, db=%g \n", da, db) ; \*/\
		U[AINDEX(i, col1, ldu)] = da ;\
		U[AINDEX(i, col2, ldu)] = db ;\
    }\
}

/*#define APPLY_GIVENS_TO_COL_M(U, ldu, col1, col2, si, ei, c, s, i, da, db, dtemp, da1, db1, dtemp1, n) { \
    ASSERT(si >=0 && si < n, "") ; \
    ASSERT(ei >=0 && ei < n, "") ; \
    i = si ; \
    for (; i < ei ; i+=2) \
    {\
        ASSIGN_TO_SCALAR(da, U, AINDEX(i, col1, ldu)) ; \
        ASSIGN_TO_SCALAR(db, U, AINDEX(i, col2, ldu)) ; \
        ASSIGN_TO_SCALAR(da1, U, AINDEX(i+1, col1, ldu)) ; \
        ASSIGN_TO_SCALAR(db1, U, AINDEX(i+1, col2, ldu)) ; \
        APPLY_GIVENS(da, db, c, s, dtemp, conjs) ; \
        APPLY_GIVENS(da1, db1, c, s, dtemp1, conjs1) ; \
        ASSIGN_TO_MATRIX(da, U, AINDEX(i, col1, ldu)) ; \
        ASSIGN_TO_MATRIX(db, U, AINDEX(i, col2, ldu)) ; \
        ASSIGN_TO_MATRIX(da1, U, AINDEX(i+1, col1, ldu)) ; \
        ASSIGN_TO_MATRIX(db1, U, AINDEX(i+1, col2, ldu)) ; \
    }\
    if (i == ei)\
    {\
        ASSIGN_TO_SCALAR(da, U, AINDEX(i, col1, ldu)) ; \
        ASSIGN_TO_SCALAR(db, U, AINDEX(i, col2, ldu)) ; \
        APPLY_GIVENS(da, db, c, s, dtemp, conjs) ; \
        ASSIGN_TO_MATRIX(da, U, AINDEX(i, col1, ldu)) ; \
        ASSIGN_TO_MATRIX(db, U, AINDEX(i, col2, ldu)) ; \
    }\
}*/

/* ========================================================================= */
/* ======= Definitions for complex/real usage  ============================= */
/* ========================================================================= */

#ifdef PIROBAND_COMPLEX /* [ */

#define MULT(c, a, b) (c[0]) = (a[0]) * (b[0]) - (a[1]) * (b[1]) ; \
		      (c[1]) = (a[0]) * (b[1]) + (a[1]) * (b[0]) ;

#define MULT_ADD(c, a, b) (c[0]) += (a[0]) * (b[0]) - (a[1]) * (b[1]) ; \
			  (c[1]) += (a[0]) * (b[1]) + (a[1]) * (b[0]) ;

#define MULT_SUB(c, a, b) (c[0]) -= ((a[0]) * (b[0]) - (a[1]) * (b[1])) ; \
			  (c[1]) -= ((a[0]) * (b[1]) + (a[1]) * (b[0])) ;

#define PRINT_NUMVALUES(str, da, db, c, s) PRINT((\
       "%s da=%0.4f %0.4fi, db=%0.4f %0.4fi, c=%0.4f%0.4fi, s=%0.4f %0.4fi\n", \
            str, da[0], da[1], db[0], db[1], c[0], c[1], s[0], s[1])) ;\

#define CSIZE(a) (a) * 2 

#define SCALE(f, scalar) f[0] = f[0] * scalar ; \
                          f[1] = f[1] * scalar ;

#define ABS(f) (((f) < 0) ? -(f) : (f)) 

#else /* ] !PIROBAND_COMPLEX  [ */

#define MULT(c, a, b) (c[0]) = (a[0]) * (b[0]) ;

#define MULT_ADD(c, a, b) (c[0]) += (a[0]) * (b[0]) ;

#define MULT_SUB(c, a, b) (c[0]) -= (a[0]) * (b[0]) ;

#define PRINT_NUMVALUES(str, da, db, c, s) PRINT((\
       "%s da=%0.4f, db=%0.4f, c=%0.4f, s=%0.4f \n", \
            str, da[0], db[0], c[0], s[0])) ;\

#define CSIZE(a) (a) 

#define SCALE(f, scalar) f[0] = f[0] * scalar ; 

#define ABS(f) (((f) < 0) ? -(f) : (f)) 

#endif /* ] PIROBAND_COMPLEX */

/* ========================================================================= */
/* === Definitions to assign scalar variables from/to mixed complex arrays = */
/* ========================================================================= */
#ifdef PIROBAND_COMPLEX /* [ */

#define ASSIGN_TO_SCALAR(s, arr, index) s[0] = arr[(index)] ; \
					s[1] = arr[(index)+1] ;

#define ASSIGN_TO_ZERO(s) s[0] = ZERO ; \
			  s[1] = ZERO ;

#define NEG(s) s[0] = -s[0] ; \
	       s[1] = -s[1] ;

#define CONJ(c, a) c[0] = a[0] ; \
                   c[1] = -a[1] ;

#define ASSIGN_TO_MATRIX(s, arr, index) arr[(index)] = s[0] ; \
					arr[(index)+1] = s[1] ;

#define ASSIGN_ZERO_TO_MATRIX(arr, index) arr[(index)] = ZERO ; \
					arr[(index)+1] = ZERO ;

#define ASSIGN_MATRIX_TO_MATRIX(s, sindex, a, index) s[(sindex)] = a[(index)] ;\
					s[(sindex)+1] = a[(index)+1] ;

#else /* ] !PIROBAND_COMPLEX  [ */

#define ASSIGN_TO_SCALAR(s, arr, index) s[0] = arr[(index)] ;

#define ASSIGN_TO_ZERO(s) s[0] = ZERO ;

#define NEG(s) s[0] = -s[0] ;

#define CONJ(c, a) c[0] = a[0] ;

#define ASSIGN_TO_MATRIX(s, arr, index) arr[(index)] = s[0] ;

#define ASSIGN_ZERO_TO_MATRIX(arr, index) arr[(index)] = ZERO ;

#define ASSIGN_MATRIX_TO_MATRIX(s, sindex, a, index) s[(sindex)] = a[(index)] ;


#endif /* ] PIROBAND_COMPLEX */

/* ============ Definitions to save the Givens rotations =================== */
#define ASSIGN_CS_SCALARS(c, s, givens, cnt) { \
		    ASSIGN_TO_SCALAR(c, givens, CSIZE(cnt*2)) ;\
		    ASSIGN_TO_SCALAR(s, givens, CSIZE(cnt*2+1)) ; }

#define SAVE_CS_SCALARS(givens, cnt, c, s) { \
		    ASSIGN_TO_MATRIX(c, givens, CSIZE(cnt*2)) ;\
		    ASSIGN_TO_MATRIX(s, givens, CSIZE(cnt*2+1)) ; }


/* ========================================================================= */
/* ======= Definitions to index a Band Matrix stored in packed form  ======  */
/* ========================================================================= */
#define INDEX(row, col) (CSIZE((row)-(col)+(obu)+((ldab)*(col)))) 

/* Definition to index (row, col) element in a dense matrix */
#define AINDEX(row, col, ld) (CSIZE((col)*(ld) + (row)))

/* ----------- Define the datatype for the update of U and V --------------- */
/* There are only two versions - for long and integer data types. */
#ifndef L_UV_HACK
typedef struct GENBAND_LONG_NAME(uv_update_struct) GENBAND_LONG_NAME(uv_update);
#define L_UV_HACK
#endif

/* ========================================================================= */
/* ==================== Definitions for internal routines  ================= */
/* ========================================================================= */
void GENBAND(reduce_equalbw)
(
    Entry *A,	            /* Band Matrix	                       */
    Int ldab,		    /* leading dimension of A                  */
    Int m,		    /* #rows in the original matrix            */
    Int n,		    /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,		    /* lower bandwidth                         */
    Int bu, 		    /* upper bandwidth                         */
    Int nc,                 /* #columns in the upper block             */
    Int nr,                 /* #rows in the upper block                */
    Int ncl,                /* #columns in the lower block             */
    Int nrl,                /* #rows in the lower block                */
    Entry *dws,             /* workspace of size 2*MAX(nc*nr, ncl*nrl) */
    Entry *B1,              /* o/p diagonal 0                          */
    Entry *B2,              /* o/p diagonal 1                          */
    Entry *U,               /* o/p accumalated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumalated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym                 /* Is the matrix symmetric ?               */
) ;

void GENBAND(reduce_unequalbw)
(
    Entry *A,	            /* Band Matrix	                       */
    Int ldab,		    /* leading dimension of A                  */
    Int m,		    /* #rows in the original matrix            */
    Int n,		    /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,		    /* lower bandwidth                         */
    Int bu, 		    /* upper bandwidth                         */
    Int nc,                 /* #columns in the upper block             */
    Int nr,                 /* #rows in the upper block                */
    Int ncl,                /* #columns in the lower block             */
    Int nrl,                /* #rows in the lower block                */
    Entry *dws,             /* workspace of size 2*MAX(nc*nr, ncl*nrl) */
    Entry *B1,              /* o/p diagonal 0                          */
    Entry *B2,              /* o/p diagonal 1                          */
    Entry *U,               /* o/p accumalated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumalated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym                 /* Is the matrix symmetric ?               */
) ;

void GENBAND(reduce_blk_lower_band)
(
    Entry *A,	            /* Band Matrix	                       */
    Int ldab,		    /* leading dimension of A                  */
    Int m,		    /* #rows in the original matrix            */
    Int n,		    /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,		    /* lower bandwidth                         */
    Int bu, 		    /* upper bandwidth                         */
    Int ncl,                /* #columns in the lower block             */
    Int nrl,                /* #rows in the lower block                */
    Entry *givens,          /* workspace for the rotations             */
    Int k,                  /* current column                          */
    Int obu,		    /* orig. lower bandwidth of A, hidden usage :A() */
    Entry *U,               /* o/p accumalated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumalated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    GENBAND_LONG_NAME(uv_update) *upd  /* Structure for U and VT update */
) ;

void GENBAND(reduce_blk_upper_band)
(
    Entry *A,	            /* Band Matrix	                       */
    Int ldab,		    /* leading dimension of A                  */
    Int m,		    /* #rows in the original matrix            */
    Int n,		    /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,		    /* lower bandwidth                         */
    Int bu, 		    /* upper bandwidth                         */
    Int nc,                 /* #columns in the upper block             */
    Int nr,                 /* #rows in the upper block                */
    Entry *givens,          /* workspace for the rotations             */
    Int k,                  /* current row                             */
    Int obu,		    /* orig. lower bandwidth of A, hidden usage :A() */
    Entry *U,               /* o/p accumalated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumalated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym,                /* Is the matrix symmetric ?               */
    GENBAND_LONG_NAME(uv_update) *upd  /* Structure for U and VT update */
) ;

void GENBAND(apply_rblk)
(
    Int m,              /* #rows in the matrix              */
    Int scol,           /* starting column of R block       */
    Int ecol,           /* end column of R block            */
    Int rgcnt,          /* #rotations to be applied         */
    Int rseed,          /* row seed for the R block         */
    Int dsrow,          /* end row for the first wave       */
    Int derow,          /* used only for ASSERT             */
    Int obu,		/* original lower b/w of A, usage hidden : A() macro */
    Int ldab,		/* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix		                   */
    Entry *givens       /* workspace for the rotations             */
) ;

void GENBAND(apply_cblk)
(
    Int n,              /* #columns in the matrix           */
    Int scol,           /* start column for the first wave  */
    Int ecol,           /* end column for the first wave    */
    Int cgcnt,          /* #rotations to be applied         */
    Int srow,           /* start column for the C block     */
    Int erow,           /* end column for the C block       */
    Int becol,          /* used only for ASSERT             */
    Int obu,		/* original lower b/w of A, usage hidden : A() macro */
    Int ldab,		/* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix		            */
    Entry *givens       /* workspace for the rotations      */
) ;

Int GENBAND(apply_fblk)
(
    Int m,              /* #rows in the matrix              */
    Int n,              /* #columns in the matrix           */
    Int srow,           /* start row for the first wave     */
    Int erow,           /* end row for the first wave       */
    Int rgcnt,          /* #rotations to be applied         */
    Int derow,          /* end row of the F block           */
    Int bscol,          /* start column of the F block      */
    Int bu,             /* upper bandwidth of the matrix    */
    Int obu,		/* original lower b/w of A, usage hidden : A() macro */
    Int ldab,		/* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix		            */
    Entry *givens,      /* workspace for the rotations      */
    Entry *V,           /* r.h.s of A to apply the rotations*/
    Int ldv,            /* leading dimension of v           */
    GENBAND_LONG_NAME(uv_update) *upd  /* Structure for U and VT update */
) ;

Int GENBAND(apply_dblk)
(
    Int m,              /* #rows in the matrix              */
    Int n,              /* #columns in the matrix           */
    Int nrc,            /* #columns in C                    */
    Int scol,           /* start column of the first wave   */
    Int ecol,           /* end column of the first wave     */
    Int cgcnt,          /* #rotations to be applied         */
    Int dsrow,          /* start row of the D block         */
    Int becol,          /* end column of the D block        */
    Int bl,             /* lower bandwidth of the matrix    */
    Int obu,		/* original lower b/w of A, usage hidden : A() macro */
    Int ldab,		/* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix		            */
    Entry *givens,      /* workspace for the rotations      */
    Entry *U,           /* o/p accumalated left rotations          */
    Int ldu,            /* leading dimension of u                  */
    Entry *C,           /* o/p for left rotations                  */
    Int ldc,            /* leading dimension of C                  */
    Int sym,            /* Is the matrix symmetric ?               */
    GENBAND_LONG_NAME(uv_update) *upd  /* Structure for U and VT update */
) ;

Int GENBAND(sblk_upper_band)
(
    Int n,              /* #columns in the matrix     */
    Int scol,		/* start column for the s-blk */
    Int ecol,		/* end column for the s-blk   */
    Int k,		/* current row                */
    Int nrow,		/* nr for the current block   */
    Int becol,		/* only for debug */
    Int serow,		/* only for debug */
    Int cseed,		/* column seed for the s-blk   */
    Int obu,		/* original lower bandwidth  of A, hidden usage :A() */
    Int ldab,		/* leading dim of the A, usage hidden in A() macro */
    Entry *A,           /* Band matrix */
    Entry *givens,      /* workspace for the rotations */
    Entry *V,           /* r.h.s of A to apply the rotations */
    Int ldv,            /* leading dimension of v */
    GENBAND_LONG_NAME(uv_update) *upd  /* Structure for U and VT update */
) ;

Int GENBAND(sblk_lower_band)
(
    Int m,		/* #rows in the matrix      */
    Int nrc,            /* #columns in C            */
    Int srow,		/* start row of s-blk       */
    Int erow,		/* end row of s-blk         */
    Int k,		/* current column           */
    Int ncol,		/* #rows in the blk         */
    Int rseed,		/* row seed for the s-blk   */
    Int obu,		/* original lower b/w of A, usage hidden : A() macro */
    Int ldab,		/* leading dim of A, usage hidden : A() macro */
    Entry *A,           /* band matrix		                   */
    Entry *givens,      /* workspace for the rotations             */
    Entry *U,           /* o/p accumalated left rotations          */
    Int ldu,            /* leading dimension of u                  */
    Entry *C,           /* o/p for left rotations                  */
    Int ldc,            /* leading dimension of C                  */
    GENBAND_LONG_NAME(uv_update) *upd  /* Structure for U and VT update */
) ;


#ifdef PIROBAND_COMPLEX

/* To find the absolute value of a complex number */
Entry GENBAND(hypot)
(
    Entry x,
    Entry y
) ;

void GENBAND(complex_to_real)
(
    Entry *A,	            /* Band Matrix	                       */
    Int ldab,		    /* leading dimension of A                  */
    Int m,		    /* #rows in the original matrix            */
    Int n,		    /* #columns in the original matrix         */
    Int nrc,                /* #columns in C                           */
    Int bl,		    /* lower bandwidth                         */
    Int bu, 		    /* upper bandwidth                         */
    Entry *U,               /* o/p accumalated left rotations          */
    Int ldu,                /* leading dimension of u                  */
    Entry *V,               /* o/p accumalated right rotations         */
    Int ldv,                /* leading dimension of v                  */
    Entry *C,               /* o/p for left rotations                  */
    Int ldc,                /* leading dimension of C                  */
    Int sym                 /* Is the matrix symmetric ?               */
) ;

#endif
