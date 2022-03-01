//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlockedKIJ.h"

void DgemmVectorNonBlockedKIJ::vector_dgemm(int n, int adjustN, int vN, const Vec4d *vA, const Vec4d *vB, Vec4d *vC) {
    for (int k = 0; k < adjustN; ++k)
        for (int i = 0; i < vN; ++i)
        {
            /* Pre-store vA[i][k] to avoid repeated memory access */
            Vec4d aik = vA[i+k*vN];

            /* Compute next partial sum of vC(i,j) */
            for( int j = 0; j < n; j++ )
                vC[i+j*vN] += aik * vB[k+j*adjustN];
        }
}

const char *DgemmVectorNonBlockedKIJ::dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with k-i-j loop order";
}
