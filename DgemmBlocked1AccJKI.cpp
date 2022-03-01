//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmBlocked1AccJKI.h"

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
void DgemmBlocked1AccJKI::do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C)
{
    for ( int j = 0; j < N; j++ )
        for (int k = 0; k < K; ++k)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bkj = A[k+j*lda];

            /* Compute next partial sum on C(i,j) */
            for(int i = 0; i < M; ++i)
                C[i+j*lda] += A[i+k*lda] * bkj;
        }
}

const char *DgemmBlocked1AccJKI::dgemm_desc() {
    return "Blocked DGEMM with 1 accumulator and j-k-i loop order inside block";
}