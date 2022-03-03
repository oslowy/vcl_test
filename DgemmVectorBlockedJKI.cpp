//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorBlockedJKI.h"

const char *dgemm_desc() {
    return "Blocked Vectorized dgemm with j-k-i loop order within block";
}

void do_block(int vM, int n, int M, int N, int K, const Vec4d *vA, const double *B, Vec4d *vC) {
    for ( int j = 0; j < N; j++ )
        for (int k = 0; k < K; ++k)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bkj = B[k + j * n];

            /* Compute next partial sum on C(i,j) */
            for(int i = 0; i < M; ++i)
                vC[i + j * vM] += vA[i + k * vM] * bkj;
        }
}
