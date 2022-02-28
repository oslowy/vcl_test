//
// Created by Serendipity_2 on 2/24/2022.
//

#include "DgemmNonBlocked1AccKIJ.h"

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void DgemmNonBlocked1AccKIJ::square_dgemm(int n, const double *A, const double *B, double *C)
{
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i)
        {
            /* Pre-store A[i][k] to avoid repeated memory access */
            double aik = A[i+k*n];

            /* Compute C(i,j) */
            for( int j = 0; j < n; j++ )
                C[i+j*n] += aik * B[k+j*n];
        }
}

const char *DgemmNonBlocked1AccKIJ::dgemm_desc() {
    return "Non-blocked DGEMM with k-i-j loop order";
}