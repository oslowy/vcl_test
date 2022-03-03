//
// Created by Orion on 2/28/2022.
//

#include "DgemmNonBlocked4AccJKI.h"

void square_dgemm(int n, const double *A, const double *B, double *C) {
    const int FOUR_ACC_LIMIT = n - 3; //Stop accumulating with four accumulators early so does not go out of bounds

    /* Compute the partial sum for four accumulators */
    for (int j = 0; j < n; ++j)
    {
        int k;
        for (k = 0; k < FOUR_ACC_LIMIT; k += 4)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bk0j = B[k+j*n];
            double bk1j = B[k+1+j*n];
            double bk2j = B[k+2+j*n];
            double bk3j = B[k+3+j*n];

            /* Compute next partial sum on C(i,j) */
            for( int i = 0; i < n; i++ )
            {
                C[i+j*n] += A[i+(k)*n] * bk0j;
                C[i+j*n] += A[i+(k+1)*n] * bk1j;
                C[i+j*n] += A[i+(k+2)*n] * bk2j;
                C[i+j*n] += A[i+(k+3)*n] * bk3j;
            }
        }

        /* Finish remaining elements not covered by the four accumulators using one accumulator */
        for(; k < n; ++k)
        {
            double bkj = B[k+j*n];

            for( int i = 0; i < n; i++ )
                C[i+j*n] += A[i+k*n] * bkj;
        }
    }
}

const char *dgemm_desc() {
    return "Non-blocked DGEMM with j-k-i loop order and four accumulators";
}
