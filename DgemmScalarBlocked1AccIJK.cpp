//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmScalarBlocked1AccIJK.h"

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
void do_block (const int n, const int M, const int N, const int K, const double* A, const double* B, double* C)
{
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Compute C(i,j) */
            double cij = C[i+ j * n]; //Need to retrieve here because each block is only a partial sum
            for (int k = 0; k < K; ++k)
                cij += A[i+ k * n] * B[k + j * n];
            C[i+ j * n] = cij;
        }
}

const char *dgemm_desc() {
    return "Blocked DGEMM with 1 accumulator and i-j-k loop order inside block";
}
