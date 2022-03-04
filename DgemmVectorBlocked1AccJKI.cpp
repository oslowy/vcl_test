//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorBlocked1AccJKI.h"

const char *dgemm_desc() {
    return "Blocked Vectorized dgemm with 1 accumulator and j-k-i loop order within block";
}

void do_block (const int n, const int M, const int N, const int K, const double* A, const double* B, double* C)
{
    const int cutM = M & (-VEC_SIZE); //Round down to multiple of 4
    const int Mremainder = M - cutM;

    for ( int j = 0; j < N; j++ )
        for (int k = 0; k < K; ++k)
        {
            /* Pre-store B(k, j) to avoid repeated memory access */
            double bkj = B[k + j * n];

            /* Compute next partial sum of vC(i,j) */
            int i;
            /* Use full load/store for the largest vector-multiple subset of the matrix */
            for(i = 0; i < cutM; i += VEC_SIZE)
                (Vec4d().load(C + i + j * n) + Vec4d().load(A + i + k * n) * bkj)
                        .store(C + i + j * n);

            /* Use partial load/store on the rest of the matrix */
            if(i < M)
                (Vec4d().load_partial(Mremainder,C + i + j * n)
                        + Vec4d().load_partial(Mremainder, A + i + k * n)
                        * bkj)
                    .store_partial(Mremainder, C + i + j * n);
        }
}
