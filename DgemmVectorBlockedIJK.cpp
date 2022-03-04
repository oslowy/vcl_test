//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorBlockedIJK.h"

const char *dgemm_desc() {
    return "Blocked Vectorized dgemm with i-j-k loop order within block";
}

void do_block (int n, int M, int N, int K, const double* A, const double* B, double* C)
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
    for (int j=0; j < N; j++)
    {
        Vec4d cij = Vec4d().load_partial(Mremainder, C + i + j * n);
        for (int k=0; k < K; k++)
            cij += Vec4d().load_partial(Mremainder,A + i + k * n) * B[k + j * n];
        cij.store_partial(Mremainder,C + i + j * n);
    }
}
