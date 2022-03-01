//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlockedJKI.h"

void DgemmVectorNonBlockedJKI::vector_dgemm(int n, int adjustN, int vN, const Vec4d *vA, const Vec4d *vB, Vec4d *vC) {
    for ( int j = 0; j < n; j++ )
        for (int k = 0; k < adjustN; ++k)
        {
            /* Pre-store vA[i][k] to avoid repeated memory access */
            Vec4d bkj = vB[k+j*adjustN];

            /* Compute next partial sum of vC(i,j) */
            for(int i = 0; i < vN; ++i)
                vC[i+j*vN] += vA[i+k*vN] * bkj;
        }
}

const char *DgemmVectorNonBlockedJKI::dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with j-k-i loop order";
}
