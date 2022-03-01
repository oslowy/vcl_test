//
// Created by Orion on 2/28/2022.
//

#include "DgemmBlocked4AccIJK.h"

const char *DgemmBlocked4AccIJK::dgemm_desc() {
    return "Blocked DGEMM with 4 accumulators and i-j-k loop order inside block";
}

void DgemmBlocked4AccIJK::do_block(int lda, int M, int N, int K, const double *A, const double *B, double *C) {
    const int FOUR_ACC_LIMIT = K - 3; //Stop accumulating with four accumulators early so does not go out of bounds

    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Initialize the accumulators */
            double acc0 = 0.0;
            double acc1 = 0.0;
            double acc2 = 0.0;
            double acc3 = 0.0;

            /* Compute the partial sums for four accumulators */
            int k;
            for(k = 0; k < FOUR_ACC_LIMIT; k += 4)
            {
                acc0 += A[i+(k)*lda] * B[k+j*lda];
                acc1 += A[i+(k+1)*lda] * B[k+1+j*lda];
                acc2 += A[i+(k+2)*lda] * B[k+2+j*lda];
                acc3 += A[i+(k+3)*lda] * B[k+3+j*lda];
            }

            /* Finish remaining elements not covered by the four accumulators using one accumulator */
            for(; k < K; k++)
                acc0 += A[i+k*lda] * B[k+j*lda];

            /* Sum the partial sums using specified independent pair association */
            C[i+j*lda] = (acc0 + acc1) + (acc2 + acc3);
        }
}
