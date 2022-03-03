//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMMVECTOR_H
#define VCL_TEST_DGEMMVECTOR_H

#include "Dgemm.h"
#include "vectorclass.h"
#include "dgemm_utils.h"

#ifndef VEC_SIZE
#define VEC_SIZE 4
#endif

void vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC);

void pad_scalar_mat(int padN, int n, double *padA, const double *A) {
    /* Copy existing data */
    for(int i=0; i < n; i++)
        for(int j=0; j < n; j++)
        {
            padA[i + j * padN] = A[i + j * n];
        }

    /* Add padding rows to A and B and padding columns to A */
    for(int i=n; i < padN; i++)
    {
        for(int j=0; j < padN; j++)
            padA[i + j * padN] = 0.0;
    }
}

void load_vA(int vM, int padN, int n, Vec4d *vA, double *padA) {
    for(int i=0; i < vM; i++) //Row of vA
        for(int j=0; j < n; j++) //Column of vA (incremented by block of VEC_SIZE)
        {
            vA[i + j * vM].load(padA + i * VEC_SIZE + j * padN);
        }
}

void load_vC(int vM, int n, Vec4d *vC) {
    for(int i=0; i<vM; i++)
        for(int j=0; j < n; j++)
            vC[i + j * vM] = 0.0;
}

void load_vectors(int vM, int n, const double *A, Vec4d *vA, Vec4d *vC) {
    /* Pad the scalar array with zeros if its size is not a multiple of the vector size */
    int padN = vM * VEC_SIZE;
    double padA[padN * n];
    /* padA: padN rows x n columns
     * padC: padN rows x n columns */

    /* Pad the matrices with zeros to evenly fill vectors */
    pad_scalar_mat(padN, n, padA, A);

    /* Load the vector version of A by expanding the cyclic permutation of each chunk of VEC_SIZE adjacent elements of A */
    load_vA(vM, padN, n, vA, padA);

    /* Initialize vC */
    load_vC(vM, n, vC);
}

void store_vectors(int vM, int n, double *C, const Vec4d *vC) {
    int padN = vM * VEC_SIZE;
    double padC[padN * n];

    /* Store padded result */
    for(int i=0; i<vM; i++)
        for(int j=0; j < n; j++)
            vC[i + j * vM].store(padC + i * VEC_SIZE + j * padN);

    /* Unpad result while copying to final result */
    for(int i=0; i < n; i++)
        for(int j=0; j < n; j++)
            C[i + j * n] = padC[i + j * padN];
}

void square_dgemm(int n, const double *A, const double *B, double *C) {
    /* Load data into vectors */
    int vM = n / VEC_SIZE + (n % VEC_SIZE == 0 ? 0 : 1); //Round up to a multiple of VEC_SIZE
    Vec4d vA[vM * n], vC[vM * n];
    /* vA: vM rows x n columns
     * vC: vM rows x n columns */

    load_vectors(vM, n, A, vA, vC);

    /* Perform blocked dgemm using vectors */
    vector_dgemm(vM, n, vA, B, vC);

    /* Extract data from vectors to final matrix */
    store_vectors(vM, n, C, vC);
}

#endif //VCL_TEST_DGEMMVECTOR_H
