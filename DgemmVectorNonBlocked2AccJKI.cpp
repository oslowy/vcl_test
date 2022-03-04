//
// Created by Orion on 3/4/2022.
//

#include "DgemmVectorNonBlocked2AccJKI.h"

/*
 * Based on example 9.12 from Agner Fog's VCL manual
 * */
void square_dgemm(const int n, const double *A, const double *B, double *C)
{
    const int TWO_ACC_LIMIT = n - 1; //Stop accumulating with two accumulators early so not out of bounds
    const int cutN = n & (-VEC_SIZE); //Round down to multiple of 4
    const int remainder = n - cutN;

    for (int j = 0; j < n; j++)
    {
        /* Compute the partial sum for four accumulators */
        int k;
        for (k = 0; k < TWO_ACC_LIMIT; k++)
        {
            /* Pre-store elements of B to avoid repeated memory access */
            double bk0j = B[k + j * n];
            double bk1j = B[k + 1 + j * n];

            /* Compute next partial sum of vC(i,j) */
            int i;
            /* Use full load/store for the largest vector-multiple subset of the matrix */
            for(i = 0; i < cutN; i += VEC_SIZE)
                (Vec4d().load(C + i + j * n) +
                    (Vec4d().load(A + i + k * n) * bk0j
                        + Vec4d().load(A + i + k + 1 * n) * bk1j))
                    .store(C + i + j * n);

            /* Use partial load/store on the rest of the matrix */
            if(i < n)
                (Vec4d().load_partial(remainder,C + i + j * n) +
                    (Vec4d().load_partial(remainder, A + i + k * n) * bk0j
                        + Vec4d().load_partial(remainder, A + i + k + 1 * n) * bk1j))
                    .store_partial(remainder, C + i + j * n);
        }

        for(; k < n; ++k)
        {
            double bkj = B[k + j * n];

            /* Compute next partial sum of vC(i,j) */
            int i;
            /* Use full load/store for the largest vector-multiple subset of the matrix */
            for(i = 0; i < cutN; i += VEC_SIZE)
                (Vec4d().load(C + i + j * n) +
                 Vec4d().load(A + i + k * n) * bkj)
                        .store(C + i + j * n);

            /* Use partial load/store on the rest of the matrix */
            if(i < n)
                (Vec4d().load_partial(remainder,C + i + j * n) +
                 Vec4d().load_partial(remainder, A + i + k * n) * bkj)
                        .store_partial(remainder, C + i + j * n);
        }
    }

}

const char *dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with 2 accumulators and j-k-i loop order";
}
