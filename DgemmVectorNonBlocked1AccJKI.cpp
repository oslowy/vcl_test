//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlocked1AccJKI.h"

void square_dgemm(const int n, const double *A, const double *B, double *C)
{
    const int cutN = n & (-VEC_SIZE); //Round down to multiple of 4
    const int remainder = n - cutN;

    for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
        {
            /* Pre-store vA[i][k] to avoid repeated memory access */
            double bkj = B[k+j*n];

            /* Compute next partial sum of vC(i,j) */
            int i;
            /* Use full load/store for the largest vector-multiple subset of the matrix */
            for(i = 0; i < cutN; i += VEC_SIZE)
                (Vec4d().load(C + i + j * n) + Vec4d().load(A + i + k * n) * bkj)
                    .store(C + i + j * n);

            /* Use partial load/store on the rest of the matrix */
            if(i < n - 1)
                (Vec4d().load_partial(remainder,C + i + j * n)
                        + Vec4d().load_partial(remainder, A + i + k * n)
                        * bkj)
                    .store_partial(remainder, C + i + j * n);
        }
}

const char *dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with j-k-i loop order";
}
