//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmNonBlocked1AccJKI.h"

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int n, const double *A, const double *B, double *C)
{
    for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bkj = B[k+j*n];

            /* Compute next partial sum on C(i,j) */
            for( int i = 0; i < n; i++ )
                C[i+j*n] += A[i+k*n] * bkj;
        }
}

const char *dgemm_desc() {
    return "Non-blocked DGEMM with j-k-i loop order and one accumulator";
}