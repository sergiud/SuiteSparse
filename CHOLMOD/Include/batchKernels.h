/* ========================================================================== */
/* === Include/batchKernels.h =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/batchKernels.h.
 * Copyright (C) 2016, Timothy A. Davis
 * CHOLMOD/Include/batchKernels.h and the CHOLMOD GPU Module are licensed under
 * Version 2.0 of the GNU General Public License.  See gpl.txt for a text of
 * the license.  CHOLMOD is also available under other licenses; contact
 * authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */
#include <cuda_runtime.h>
#include <cublas_v2.h>

#ifdef __cplusplus
extern "C" {
#endif


/* batched dgemm */
void dgemm_custom_simple_1block_batch (cudaStream_t stream,
                                       cublasOperation_t transa,
                                       cublasOperation_t transb,
                                       int *mlist, int *nlist, int *klist,
                                       const double *alpha,
                                       const double **Alist, int *ldalist,
                                       const double **Blist, int *ldblist,
                                       const double *beta,
                                       double **Clist, int *ldclist, int nbatch);

/* batched dsyrk */
void dsyrk_custom_simple_1block_batch(cudaStream_t stream,
                                      cublasFillMode_t uplo,
                                      cublasOperation_t transa,
                                      int *nlist, int *klist,
                                      const double *alpha,
                                      const double **Alist, int *ldalist,
                                      const double *beta,
                                      double **Clist, int *ldclist, int nbatch);

/* batched dtrsm */
void dtrsm_custom_simple_1block_batch (cudaStream_t stream,
                                       cublasSideMode_t side,
                                       cublasFillMode_t uplo,
                                       cublasOperation_t trans,
                                       cublasDiagType_t diag,
                                       int *mlist, int *nlist,
                                       const double *alpha,
                                       const double **Alist, int *ldalist,
                                       double **Blist, int *ldblist, int nbatch);

/* batched dpotrf */
void dpotrf_custom_simple_1block_batch(cudaStream_t stream,
                                       cublasFillMode_t uplo,
                                       int *nlist,
                                       double **alist,
                                       int *ldalist, int *info, int nbatch);


/* custom CPU potrf functions */
int dpotf2_custom_cpu(char *uplo, int *n, double *a, int *lda, int *info);
int dpotrf_custom_cpu(char *uplo, int *n, double *a, int *lda, int *info);


#ifdef __cplusplus
}
#endif

