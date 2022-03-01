//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlockedIJK.h"

void DgemmVectorNonBlockedIJK::vector_dgemm(int vM, int vN, int vK, const Vec4d *vA, const Vec4d *vB, Vec4d *vC) {
    /* For each row i of vA */
    for (int i = 0; i < vM; ++i)
        /* For each column j of vB */
        for (int j = 0; j < vK; ++j)
        {
            /* Compute vC(i,j) */
            Vec4d cij = vC[i+j*vM];
            for( int k = 0; k < vN; k++ )
                cij += vA[i+k*vM] * vB[k+j*vN];
            vC[i+j*vM] = cij;
        }
}

const char *DgemmVectorNonBlockedIJK::dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with i-j-k loop order";
}
