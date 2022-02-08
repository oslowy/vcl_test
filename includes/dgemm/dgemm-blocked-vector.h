//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMM_BLOCKED_H
#define VCL_TEST_DGEMM_BLOCKED_H
#endif

#include "vectorclass.h"

const char* dgemm_desc = "Vectorized blocked dgemm.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 32
#endif

#ifndef VEC_SIZE
#define VEC_SIZE 4

#define dgemm_blocked_min(a,b) (((a)<(b))?(a):(b))

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{
    /* Initialize lookup index vectors */
    Vec4q iInit;
    for(int i = 0; i < VEC_SIZE; i++) {
        iInit.insert(i, i * lda);
    }

    Vec4q iA(iInit);
    Vec4q iB(iInit);
    Vec4q iC(iInit);

    /*
     * Need to work on remainders that do not fill size-4 vectors!
     * */

    /* For each row i of A */
    for (int i = 0; i < M; ++i) {
        /* For each column j of B */
        for (int j = 0; j < N; j += VEC_SIZE)
        {
            /* Compute C(i,j) */
//            double cij = C[i+j*lda];
            Vec4d vC = lookup<VEC_SIZE>(iC, C);

            for (int k = 0; k < K; k += VEC_SIZE)
//                cij += A[i+k*lda] * B[k+j*lda];
            {
                Vec4d vA = lookup<VEC_SIZE>(iA, A);
                Vec4d vB = lookup<VEC_SIZE>(iB, B);
                vC += vA * vB;
                vC.store(C + iC.extract(0));

                iA += VEC_SIZE * lda;
                iB += VEC_SIZE;
            }
            iB += VEC_SIZE * lda - K;
            iC += VEC_SIZE * lda;
        }
        iA += VEC_SIZE - lda * K;
        iB = iInit;
        iC += VEC_SIZE;
    }



}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm (int lda, double* A, double* B, double* C)
{
    /* For each block-row of A */
    for (int i = 0; i < lda; i += BLOCK_SIZE)
        /* For each block-column of B */
        for (int j = 0; j < lda; j += BLOCK_SIZE)
            /* Accumulate block dgemms into block of C */
            for (int k = 0; k < lda; k += BLOCK_SIZE)
            {
                /* Correct block dimensions if block "goes off edge of" the matrix */
                int M = dgemm_blocked_min (BLOCK_SIZE, lda - i);
                int N = dgemm_blocked_min (BLOCK_SIZE, lda - j);
                int K = dgemm_blocked_min (BLOCK_SIZE, lda - k);

                /* Perform individual block dgemm */
                do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
            }
}

#endif //VCL_TEST_DGEMM_BLOCKED_H
