//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmScalarBlocked1AccJKI.h"

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
void do_block (const int n, const int M, const int N, const int K, const double* A, const double* B, double* C)
{
    for ( int j = 0; j < N; j++ )
        for (int k = 0; k < K; ++k)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bkj = B[k+ j * n];

            /* Compute next partial sum on C(i,j) */
            for(int i = 0; i < M; ++i)
                C[i+ j * n] += A[i + k * n] * bkj;
        }
}

const char *dgemm_desc() {
    return "Blocked DGEMM with 1 accumulator and j-k-i loop order inside block";
}