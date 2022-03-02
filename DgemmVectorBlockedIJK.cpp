//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorBlockedIJK.h"

const char *DgemmVectorBlockedIJK::dgemm_desc() {
    return "Blocked Vectorized dgemm with i-j-k loop order within block";
}

void DgemmVectorBlockedIJK::do_block(int vM, int n, int M, int N, int K, const Vec4d *vA, const double *B, Vec4d *vC) {
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Compute C(i,j) */
            Vec4d cij = vC[i + j * vM];
            for (int k = 0; k < K; ++k)
                cij += vA[i + k * vM] * B[k + j * n];
            vC[i + j * vM] = cij;
        }
}
