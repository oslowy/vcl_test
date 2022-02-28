//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmBlocked1AccKIJ.h"

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
void DgemmBlocked1AccKIJ::do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C)
{
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < M; ++i)
        {
            /* Pre-store A[i][k] to avoid repeated memory access */
            double aik = A[i+k*lda];

            /* Compute C(i,j) */
            for( int j = 0; j < N; j++ )
                C[i+j*lda] += aik * B[k+j*lda];
        }
}

const char *DgemmBlocked1AccKIJ::dgemm_desc() {
    return "Blocked DGEMM with 1 accumulator and k-i-j loop order inside block";
}