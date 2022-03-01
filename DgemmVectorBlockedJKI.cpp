//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorBlockedJKI.h"

const char *DgemmVectorBlockedJKI::dgemm_desc() {
    return "Blocked Vectorized dgemm with j-k-i loop order within block";
}

void DgemmVectorBlockedJKI::do_block(int vM, int vN, int vK, int M, int N, int K, const Vec4d *vA, const Vec4d *vB,
                                     Vec4d *vC) {
    for ( int j = 0; j < N; j++ )
        for (int k = 0; k < K; ++k)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            Vec4d bkj = vB[k + j * vN];

            /* Compute next partial sum on C(i,j) */
            for(int i = 0; i < M; ++i)
                vC[i + j * vM] += vA[i + k * vM] * bkj;
        }
}
