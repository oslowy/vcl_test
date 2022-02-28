//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlockedIJK.h"

void DgemmVectorNonBlockedIJK::vector_dgemm(int n, int adjustN, int vN, const Vec4d *vA, const Vec4d *vB, Vec4d *vC) {
    /* For each row i of vA */
    for (int i = 0; i < n; ++i)
        /* For each column j of vB */
        for (int j = 0; j < vN; ++j)
        {
            /* Compute vC(i,j) */
            Vec4d cij = vC[i+j*n];
            for( int k = 0; k < n; k++ )
                cij += vA[i+k*adjustN] * vB[k+j*adjustN];
            vC[i+j*n] = cij;
        }
}
