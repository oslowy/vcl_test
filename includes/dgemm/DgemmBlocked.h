//
// Created by Orion on 2/22/2022.
//

#ifndef VCL_TEST_DGEMMBLOCKED_H
#define VCL_TEST_DGEMMBLOCKED_H

#include "Dgemm.h"
#include "dgemm_utils.h"

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 36
#endif

void do_block (int n, int M, int N, int K, const double* A, const double* B, double* C);

void square_dgemm(const int n, const double *A, const double *B, double *C)
{
    /* For each block-row of A */
    for (int i = 0; i < n; i += BLOCK_SIZE)
        /* For each block-column of B */
        for (int j = 0; j < n; j += BLOCK_SIZE)
            /* Accumulate block dgemms into block of C */
            for (int k = 0; k < n; k += BLOCK_SIZE)
            {
                /* Correct block dimensions if block "goes off edge of" the matrix */
                int M = min2 (BLOCK_SIZE, n - i);
                int N = min2 (BLOCK_SIZE, n - j);
                int K = min2 (BLOCK_SIZE, n - k);

                /* Perform individual block dgemm */
                do_block(n, M, N, K, A + i + k * n, B + k + j * n, C + i + j * n);
            }
}

#endif //VCL_TEST_DGEMMBLOCKED_H
