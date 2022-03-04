//
// Created by Orion on 3/4/2022.
//

#include "DgemmVectorNonBlocked2AccIJK.h"

/*
 * Based on example 9.12 from Agner Fog's VCL manual
 * */
void square_dgemm(const int n, const double *A, const double *B, double *C)
{
    const int TWO_ACC_LIMIT = n - 1; //Stop accumulating with two accumulators early so not out of bounds
    const int cutN = n & (-VEC_SIZE); //Round down to multiple of 4
    const int remainder = n - cutN;

    /* Use full load/store for the largest vector-multiple subset of the matrix */
    int i; //Keep this value after last loop iteration for an extra ending iteration
    for(i=0; i < cutN; i += VEC_SIZE)
        for (int j=0; j < n; j++)
        {
            /* Initialize the accumulators */
            Vec4d acc0(0.0);
            Vec4d acc1(0.0);

            /* Compute the partial sum for two accumulators*/
            int k;
            for (k=0; k < TWO_ACC_LIMIT; k += 2)
            {
                acc0 += Vec4d().load(A + i + (k) * n) * B[k + j * n];
                acc1 += Vec4d().load(A + i + (k + 1) * n) * B[k + 1 + j * n];
            }

            /* Finish remaining elements not covered by the two accumulators using one accumulator */
            for(; k < n; k++)
                acc0 += Vec4d().load(A + i + (k) * n) * B[k + j * n];

            /* Add the partial sums and store */
            (acc0 + acc1).store(C + i + j * n);
        }

    /* Use partial load/store on the rest of the matrix */
    if(i < n - 1)
        for (int j=0; j < n; j++)
        {
            Vec4d acc0(0.0);
            Vec4d acc1(0.0);

            /* Compute the partial sum for two accumulators*/
            int k;
            for (k=0; k < TWO_ACC_LIMIT; k += 2)
            {
                acc0 += Vec4d().load_partial(remainder, A + i + (k) * n) * B[k + j * n];
                acc1 += Vec4d().load_partial(remainder, A + i + (k + 1) * n) * B[k + 1 + j * n];
            }

            /* Finish remaining elements not covered by the two accumulators using one accumulator */
            for(; k < n; k++)
                acc0 += Vec4d().load_partial(remainder, A + i + k * n) * B[k + j * n];

            /* Add the partial sums and store */
            (acc0 + acc1).store_partial(remainder, C + i + j * n);
        }
}

const char *dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with 2 accumulators and i-j-k loop order";
}
