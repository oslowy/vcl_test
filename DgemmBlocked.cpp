//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmBlocked.h"
#include "dgemm_utils.h"

void DgemmBlocked::square_dgemm(int lda, const double *A, const double *B, double *C) {
    /* For each block-row of A */
    for (int i = 0; i < lda; i += BLOCK_SIZE)
        /* For each block-column of B */
        for (int j = 0; j < lda; j += BLOCK_SIZE)
            /* Accumulate block dgemms into block of C */
            for (int k = 0; k < lda; k += BLOCK_SIZE)
            {
                /* Correct block dimensions if block "goes off edge of" the matrix */
                int M = min2 (BLOCK_SIZE, lda - i);
                int N = min2 (BLOCK_SIZE, lda - j);
                int K = min2 (BLOCK_SIZE, lda - k);

                /* Perform individual block dgemm */
                do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
            }
}
