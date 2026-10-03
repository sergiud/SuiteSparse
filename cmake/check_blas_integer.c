/* SPDX-FileCopyrightText: 2026 Sergiu Deitsch
   SPDX-License-Identifier: Apache-2.0 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

extern void SUITESPARSE_DCOPY(const void *n, const double *x,
    const void *incx, double *y, const void *incy);
extern void SUITESPARSE_DGETRF(const void *m, const void *n,
    double *a, const void *lda, void *ipiv, void *info);

int main(void)
{
    const int64_t one = 1;
    const int64_t negative_with_positive_low_word =
        -(int64_t)UINT32_MAX;
    const int64_t zero = 0;
    uint32_t first_word;
    const uint32_t increment_words[] = {1, 0};
    int64_t increment;
    const double source = 1.0;
    double first = 0.0;
    double second = 0.0;
    size_t blas_size;
    size_t lapack_size;
    const uint32_t leading_dimension_words[] = {1, 1};
    int64_t leading_dimension;
    int64_t pivot = 0;
    uint64_t info = UINT64_MAX;
    uint32_t info_words[2];
    double matrix = 0.0;

    memcpy(&first_word, &one, sizeof(first_word));
    memcpy(&increment, increment_words, sizeof(increment));
    /* Both calls copy at most one element for either integer ABI. On little
       endian systems, the negative 64-bit count reads as a positive 32-bit
       count. On big endian systems, the positive 64-bit count reads as zero
       for a 32-bit interface. */
    SUITESPARSE_DCOPY(&one, &source, &increment, &first, &increment);
    SUITESPARSE_DCOPY(&negative_with_positive_low_word, &source,
        &increment, &second, &increment);
    if (first == source && fpclassify(second) == FP_ZERO) {
        blas_size = sizeof(int64_t);
    } else if ((first_word == 1 && first == source && second == source) ||
               (first_word == 0 && fpclassify(first) == FP_ZERO &&
                fpclassify(second) == FP_ZERO)) {
        blas_size = sizeof(int32_t);
    } else {
        fprintf(stderr, "BLAS integer ABI probe returned actual values %g and %g, "
            "expected 1 and 0, 1 and 1, or 0 and 0\n", first, second);
        return 1;
    }

    /* A zero-sized factorization writes only INFO. The leading dimension is
       positive for both integer ABIs on either byte order. */
    memcpy(&leading_dimension, leading_dimension_words, sizeof(leading_dimension));
    SUITESPARSE_DGETRF(&zero, &zero, &matrix, &leading_dimension, &pivot, &info);
    memcpy(info_words, &info, sizeof(info_words));
    if (info == 0) {
        lapack_size = sizeof(int64_t);
    } else if (info_words[0] == 0 && info_words[1] == UINT32_MAX) {
        lapack_size = sizeof(int32_t);
    } else {
        fprintf(stderr, "LAPACK integer ABI probe returned actual words %u and %u, "
            "expected 0 and 0 or 0 and UINT32_MAX\n",
            (unsigned int)info_words[0], (unsigned int)info_words[1]);
        return 1;
    }

    if (blas_size != lapack_size) {
        fprintf(stderr, "BLAS integer size %zu differs from LAPACK integer size %zu, "
            "expected matching integer sizes\n", blas_size, lapack_size);
        return 1;
    }
    printf("%zu", blas_size);
    return 0;
}
