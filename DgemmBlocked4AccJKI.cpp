//
// Created by Orion on 2/28/2022.
//

#include "DgemmBlocked4AccJKI.h"

const char *dgemm_desc() {
    return "Blocked DGEMM with 4 accumulators and j-k-i loop order inside block";
}

void do_block(int lda, int M, int N, int K, const double *A, const double *B, double *C) {
    const int FOUR_ACC_LIMIT = K - 3; //Stop accumulating with four accumulators early so does not go out of bounds

    for ( int j = 0; j < N; j++ )
    {
        int k;
        for (k = 0; k < FOUR_ACC_LIMIT; k += 4)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bk0j = B[k+j*lda];
            double bk1j = B[k+1+j*lda];
            double bk2j = B[k+2+j*lda];
            double bk3j = B[k+3+j*lda];

            /* Compute next partial sum on C(i,j) */
            for(int i = 0; i < M; ++i)
            {
                C[i+j*lda] += A[i+(k)*lda] * bk0j;
                C[i+j*lda] += A[i+(k+1)*lda] * bk1j;
                C[i+j*lda] += A[i+(k+2)*lda] * bk2j;
                C[i+j*lda] += A[i+(k+3)*lda] * bk3j;
            }

        }

        /* Finish remaining elements not covered by the four accumulators using one accumulator */
        for(; k < K; ++k)
        {
            double bkj = B[k+j*lda];

            for( int i = 0; i < M; i++ )
                C[i+j*lda] += A[i+k*lda] * bkj;
        }
    }

}
