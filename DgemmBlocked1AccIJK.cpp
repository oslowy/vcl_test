//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmBlocked1AccIJK.h"

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
void do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C)
{
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Compute C(i,j) */
            double cij = C[i+j*lda]; //Need to retrieve here because each block is only a partial sum
            for (int k = 0; k < K; ++k)
                cij += A[i+k*lda] * B[k+j*lda];
            C[i+j*lda] = cij;
        }
}

const char *dgemm_desc() {
    return "Blocked DGEMM with 1 accumulator and i-j-k loop order inside block";
}
