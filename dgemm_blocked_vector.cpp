//
// Created by Orion on 2/22/2022.
//

#include "dgemm_blocked_vector.h"
#include "vectorclass.h"
#include "dgemm_utils.h"

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
void DgemmBlockedVector::do_block(int lda, int M, int N, int K, double *A, double *B, double *C)
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
     * TODO Need to work on remainders that do not fill size-4 vectors!
     * */

    /* For each row i of A */
    for (int i = 0; i < M; ++i) {
        /* For each column j of B */
        for (int j = 0; j < N; j += VEC_SIZE)
        {
            /* Compute C(i,j) */
//            double cij = C[i+j*lda];
            Vec4d vC = lookup<VEC_SIZE>(iC, C);

            // TODO fill partial vectors at end of row/column

            for (int k = 0; k < K; k += VEC_SIZE)
//                cij += A[i+k*lda] * B[k+j*lda];
            {
                Vec4d vA = lookup<VEC_SIZE>(iA, A);
                Vec4d vB = lookup<VEC_SIZE>(iB, B);
                vC += vA * vB;

                // TODO store only the part of the vector that was filled with real data when it is partial
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
