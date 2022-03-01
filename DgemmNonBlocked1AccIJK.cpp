//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmNonBlocked1AccIJK.h"

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void DgemmNonBlocked1AccIJK::square_dgemm(int n, const double *A, const double *B, double *C)
{
    /* For each row i of A */
    for (int i = 0; i < n; ++i)
        /* For each column j of B */
        for (int j = 0; j < n; ++j)
        {
            /* Compute C(i,j) */
            double cij = 0.0;
            for( int k = 0; k < n; k++ )
                cij += A[i+k*n] * B[k+j*n];
            C[i+j*n] = cij;
        }
}

const char *DgemmNonBlocked1AccIJK::dgemm_desc() {
    return "Non-blocked DGEMM with i-j-k loop order and one accumulator";
}