//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlockedIJK.h"

void DgemmVectorNonBlockedIJK::vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) {
    /* For each row i of vA */
    for (int i = 0; i < vM; ++i)
        /* For each column j of B */
        for (int j = 0; j < n; ++j)
        {
            /* Compute vC(i,j) */
            Vec4d cij = 0.0;
            for( int k = 0; k < n; k++ )
                cij += vA[i+k*vM] * B[k+j*n]; //Uniformly multiply vector A by scalar B
            vC[i+j*vM] = cij;
        }
}

const char *DgemmVectorNonBlockedIJK::dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with i-j-k loop order";
}
