//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorBlocked1AccIJK.h"

const char *dgemm_desc() {
    return "Blocked Vectorized dgemm with 1 accumulator and i-j-k loop order within block";
}

void do_block (const int n, const int M, const int N, const int K, const double* A, const double* B, double* C)
{
    const int cutM = M & (-VEC_SIZE); //Round down to multiple of 4
    const int Mremainder = M - cutM;

    /* Use full load/store for the largest vector-multiple subset of the matrix */
    int i; //Keep this value after last loop iteration for an extra ending iteration
    for(i=0; i < cutM; i += VEC_SIZE)
        for (int j=0; j < N; j++)
        {
            Vec4d cij = Vec4d().load(C + i + j * n);
            for (int k=0; k < K; k++)
                cij += Vec4d().load(A + i + k * n) * B[k + j * n];
            cij.store(C + i + j * n);
        }

    /* Use partial load/store on the rest of the matrix */
    if(i < M)
        for (int j=0; j < N; j++)
        {
            Vec4d cij = Vec4d().load_partial(Mremainder, C + i + j * n);
            for (int k=0; k < K; k++)
                cij += Vec4d().load_partial(Mremainder,A + i + k * n) * B[k + j * n];
            cij.store_partial(Mremainder,C + i + j * n);
        }
}
