//
// Created by Orion on 2/28/2022.
//

#include "DgemmVectorNonBlocked1AccIJK.h"

void vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) {
    /* For each row i of vA */
    for (int i = 0; i < vM; ++i)
        /* For each column j of B */
        for (int j = 0; j < n; ++j)
        {
            /* Compute vC(i,j) */
            Vec4d cij = 0.0;
            for( int k = 0; k < n; k++ )
                cij += vA[i+k*vM] * B[k+j*n]; //Uniformly multiply vector A by scalar B
            vC[i+j*vM] = cij;
        }
}

/*
 * Based on example 9.12 from Agner Fog's VCL manual
 * */
void square_dgemm(int n, const double *A, const double *B, double *C)
{
    const int cutN = n & (-VEC_SIZE); //Round down to multiple of 4
    const int remainder = n - cutN;

    /* Use full load/store for the largest vector-multiple subset of the matrix */
    int i; //Keep this value after last loop iteration for an extra ending iteration
    for(i=0; i < cutN; i += VEC_SIZE)
        for (int j=0; j < n; j++)
        {
            Vec4d cij(0.0);
            for (int k=0; k < n; k++)
                cij += Vec4d().load(A + i + k * n) * B[k + j * n];
            cij.store(C + i + j * n);
        }

    /* Use partial load/store on the rest of the matrix */
    for (int j=0; j < n; j++)
    {
        Vec4d cij(0.0);
        for (int k=0; k < n; k++)
            cij += Vec4d().load_partial(remainder,A + i + k * n) * B[k + j * n];
        cij.store_partial(remainder,C + i + j * n);
    }
}

const char *dgemm_desc() {
    return "Non-blocked Vectorized DGEMM with i-j-k loop order";
}
