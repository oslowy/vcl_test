//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorBlockedIJK.h"

const char *DgemmVectorBlockedIJK::dgemm_desc() {
    return "Blocked Vectorized dgemm with i-j-k loop order within block";
}

void DgemmVectorBlockedIJK::do_block(int vM, int vN, int vK, int M, int N, int K, const Vec4d *vA, const Vec4d *vB,
                                     Vec4d *vC) {
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Compute C(i,j) */
            Vec4d cij = 0.0;
            for (int k = 0; k < K; ++k)
                cij += vA[i + k * vM] * vB[k + j * vN];
            vC[i + j * vM] = cij;
        }
}
