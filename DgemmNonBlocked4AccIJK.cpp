//
// Created by Orion on 2/28/2022.
//

#include "DgemmNonBlocked4AccIJK.h"

void square_dgemm(int n, const double *A, const double *B, double *C) {
    const int FOUR_ACC_LIMIT = n - 3; //Stop accumulating with four accumulators early so does not go out of bounds

    /* For each row i of A */
    for (int i = 0; i < n; i++)
        /* For each column j of B */
        for (int j = 0; j < n; ++j)
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
                acc0 += A[i+(k)*n] * B[k+j*n];
                acc1 += A[i+(k+1)*n] * B[k+1+j*n];
                acc2 += A[i+(k+2)*n] * B[k+2+j*n];
                acc3 += A[i+(k+3)*n] * B[k+3+j*n];
            }

            /* Finish remaining elements not covered by the four accumulators using one accumulator */
            for(; k < n; k++)
                acc0 += A[i+k*n] * B[k+j*n];

            /* Sum the partial sums using specified independent pair association */
            C[i+j*n] = (acc0 + acc1) + (acc2 + acc3);
        }
}

const char *dgemm_desc() {
    return "Non-blocked DGEMM with i-j-k loop order and four accumulators";
}
