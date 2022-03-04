//
// Created by Orion on 2/28/2022.
//

#include "DgemmScalarBlocked4AccJKI.h"

const char *dgemm_desc() {
    return "Blocked DGEMM with 4 accumulators and j-k-i loop order inside block";
}

void do_block(const int n, const int M, const int N, const int K, const double *A, const double *B, double *C) {
    const int FOUR_ACC_LIMIT = K - 3; //Stop accumulating with four accumulators early so does not go out of bounds

    for ( int j = 0; j < N; j++ )
    {
        int k;
        for (k = 0; k < FOUR_ACC_LIMIT; k += 4)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bk0j = B[k+ j * n];
            double bk1j = B[k+1+ j * n];
            double bk2j = B[k+2+ j * n];
            double bk3j = B[k+3+ j * n];

            /* Compute next partial sum on C(i,j) */
            for(int i = 0; i < M; ++i)
            {
                C[i+ j * n] += A[i + (k) * n] * bk0j;
                C[i+ j * n] += A[i + (k + 1) * n] * bk1j;
                C[i+ j * n] += A[i + (k + 2) * n] * bk2j;
                C[i+ j * n] += A[i + (k + 3) * n] * bk3j;
            }

        }

        /* Finish remaining elements not covered by the four accumulators using one accumulator */
        for(; k < K; ++k)
        {
            double bkj = B[k+ j * n];

            for( int i = 0; i < M; i++ )
                C[i+ j * n] += A[i + k * n] * bkj;
        }
    }

}
