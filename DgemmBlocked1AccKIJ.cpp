//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmBlocked1AccKIJ.h"

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
void DgemmBlocked1AccKIJ::do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C)
{
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < M; ++i)
        {
            /* Pre-store A[i][k] to avoid repeated memory access */
            double aik = A[i+k*lda];

            /* Compute C(i,j) */
            for( int j = 0; j < N; j++ )
                C[i+j*lda] += aik * B[k+j*lda];
        }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void DgemmBlocked1AccKIJ::square_dgemm (int lda, const double* A, const double* B, double* C)
{
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

const char *DgemmBlocked1AccKIJ::dgemm_desc() {
    return "Blocked DGEMM with 1 accumulator and k-i-j loop order inside block";
}