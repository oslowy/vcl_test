//
// Created by Orion on 2/28/2022.
//

#include "DgemmNonBlocked4AccJKI.h"

void DgemmNonBlocked4AccJKI::square_dgemm(int n, const double *A, const double *B, double *C) {
    const int FOUR_ACC_LIMIT = n - 3; //Stop accumulating with four accumulators early so does not go out of bounds

    for (int j = 0; j < n; ++j)
    {
        int k;
        for (k = 0; k < FOUR_ACC_LIMIT; ++k)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bkj = B[k+j*n];

            /* Compute next partial sum on C(i,j) */
            for( int i = 0; i < n; i++ )
                C[i+j*n] += A[i+k*n] * bkj;
        }

        for(; k < n; ++k)
        {
            double bkj = A[k+j*n];

            /* Compute next partial sum on C(i,j) */
            for( int i = 0; i < n; i++ )
                C[i+j*n] += A[i+k*n] * bkj;
        }
    }

}

const char *DgemmNonBlocked4AccJKI::dgemm_desc() {
    return "Non-blocked DGEMM with j-k-i loop order and four accumulators";
}
