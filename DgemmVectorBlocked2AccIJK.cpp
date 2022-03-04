//
// Created by Orion on 3/4/2022.
//

#include "DgemmVectorBlocked2AccIJK.h"

const char *dgemm_desc() {
    return "Blocked Vectorized dgemm with 2 accumulators and i-j-k loop order within block";
}

void do_block (const int n, const int M, const int N, const int K, const double* A, const double* B, double* C)
{
    const int TWO_ACC_LIMIT = K - 1; //Stop accumulating with two accumulators early so not out of bounds
    const int cutM = M & (-VEC_SIZE); //Round down to multiple of 4
    const int Mremainder = M - cutM;

    /* Use full load/store for the largest vector-multiple subset of the matrix */
    int i; //Keep this value after last loop iteration for an extra ending iteration
    for(i=0; i < cutM; i += VEC_SIZE)
        for (int j=0; j < N; j++)
        {
            /* Initialize the accumulators */
            Vec4d acc0 = Vec4d().load(C + i + j * n);
            Vec4d acc1(0.0);

            /* Compute the partial sum for two accumulators*/
            int k;
            for (k=0; k < TWO_ACC_LIMIT; k += 2)
            {
                acc0 += Vec4d().load(A + i + (k) * n) * B[k + j * n];
                acc1 += Vec4d().load(A + i + (k + 1) * n) * B[k + 1 + j * n];
            }

            /* Finish remaining elements not covered by the two accumulators using one accumulator */
            for(; k < K; k++)
                acc0 += Vec4d().load(A + i + (k) * n) * B[k + j * n];

            /* Add the partial sums and store */
            (acc0 + acc1).store(C + i + j * n);
        }

    /* Use partial load/store on the rest of the matrix */
    if(i < M)
        for (int j=0; j < N; j++)
        {
            Vec4d acc0 = Vec4d().load_partial(Mremainder, C + i + j * n);
            Vec4d acc1(0.0);

            /* Compute the partial sum for two accumulators*/
            int k;
            for (k=0; k < TWO_ACC_LIMIT; k += 2)
            {
                acc0 += Vec4d().load_partial(Mremainder, A + i + (k) * n) * B[k + j * n];
                acc1 += Vec4d().load_partial(Mremainder, A + i + (k + 1) * n) * B[k + 1 + j * n];
            }

            /* Finish remaining elements not covered by the two accumulators using one accumulator */
            for(; k < K; k++)
                acc0 += Vec4d().load_partial(Mremainder, A + i + k * n) * B[k + j * n];

            /* Add the partial sums and store */
            (acc0 + acc1).store_partial(Mremainder, C + i + j * n);
        }
}