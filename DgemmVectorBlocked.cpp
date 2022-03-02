//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmVectorBlocked.h"
#include "dgemm_utils.h"

void DgemmVectorBlocked::vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) {
    /* For each block-row of A */
    for (int i = 0; i < vM; i += BLOCK_SIZE)
        /* For each block-column of B */
        for (int j = 0; j < n; j += BLOCK_SIZE)
            /* Accumulate block dgemms into block of C */
            for (int k = 0; k < n; k += BLOCK_SIZE)
            {
                /* Correct block dimensions if block "goes off edge of" the matrix */
                int M = min2 (BLOCK_SIZE, vM - i);
                int N = min2 (BLOCK_SIZE, n - j);
                int K = min2 (BLOCK_SIZE, n - k);

                /* Perform individual block dgemm */
                do_block(vM, n, M, N, K, vA + i + k * vM, B + k + j * n, vC + i + j * vM);
            }
}
