//
// Created by Orion on 3/4/2022.
//

#include "DgemmVectorBlocked2AccJKI.h"

const char *dgemm_desc() {
    return "Blocked Vectorized dgemm with 2 accumulators and j-k-i loop order within block";
}

void do_block (const int n, const int M, const int N, const int K, const double* A, const double* B, double* C)
{
    const int TWO_ACC_LIMIT = K - 1; //Stop accumulating with two accumulators early so not out of bounds
    const int cutM = M & (-VEC_SIZE); //Round down to multiple of 4
    const int Mremainder = M - cutM;

    for (int j = 0; j < N; j++)
    {
        /* Compute the partial sum for four accumulators */
        int k;
        for (k = 0; k < TWO_ACC_LIMIT; k += 2)
        {
            /* Pre-store elements of B to avoid repeated memory access */
            double bk0j = B[k + j * n];
            double bk1j = B[k + 1 + j * n];

            /* Compute next partial sum of C(i,j) */
            int i;
            /* Use full load/store for the largest vector-multiple subset of the matrix */
            for(i = 0; i < cutM; i += VEC_SIZE)
                (Vec4d().load(C + i + j * n) +
                 (Vec4d().load(A + i + (k) * n) * bk0j
                  + Vec4d().load(A + i + (k + 1) * n) * bk1j))
                        .store(C + i + j * n);

            /* Use partial load/store on the rest of the matrix */
            if(i < M)
                (Vec4d().load_partial(Mremainder, C + i + j * n) +
                 (Vec4d().load_partial(Mremainder, A + i + (k) * n) * bk0j
                  + Vec4d().load_partial(Mremainder, A + i + (k + 1) * n) * bk1j))
                        .store_partial(Mremainder, C + i + j * n);
        }

        for(; k < K; ++k)
        {
            double bkj = B[k + j * n];

            /* Compute next partial sum of C(i,j) */
            int i;
            /* Use full load/store for the largest vector-multiple subset of the matrix */
            for(i = 0; i < cutM; i += VEC_SIZE)
                (Vec4d().load(C + i + j * n) +
                 Vec4d().load(A + i + k * n) * bkj)
                        .store(C + i + j * n);

            /* Use partial load/store on the rest of the matrix */
            if(i < M)
                (Vec4d().load_partial(Mremainder, C + i + j * n) +
                 Vec4d().load_partial(Mremainder, A + i + k * n) * bkj)
                        .store_partial(Mremainder, C + i + j * n);
        }
    }

}