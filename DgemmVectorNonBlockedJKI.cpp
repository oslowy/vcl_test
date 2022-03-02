//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlockedJKI.h"

void DgemmVectorNonBlockedJKI::vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) {
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

const char *DgemmVectorNonBlockedJKI::dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with j-k-i loop order";
}
