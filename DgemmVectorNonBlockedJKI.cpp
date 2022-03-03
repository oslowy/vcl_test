//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlockedJKI.h"

void vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) {
    for (int j = 0; j < n; j++ )
        for (int k = 0; k < n; ++k)
        {
            /* Pre-store vA[i][k] to avoid repeated memory access */
            double bkj = B[k+j*n];

            /* Compute next partial sum of vC(i,j) */
            for(int i = 0; i < vM; ++i)
                vC[i+j*vM] += vA[i+k*vM] * bkj;
        }
}

void square_dgemm(int n, const double *A, const double *B, double *C)
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
                (Vec4d().load(A + i + k * n) * bkj).store(C + i + j * n);

            /* Use partial load/store on the rest of the matrix */
            (Vec4d().load_partial(remainder, A + i + k * n) * bkj)
                .store_partial(remainder, C + i + j * n);
        }
}

const char *dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with j-k-i loop order";
}
