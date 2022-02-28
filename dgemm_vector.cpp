//
// Created by Orion on 2/22/2022.
//

#include "dgemm_vector.h"
#include "dgemm_utils.h"

void DgemmVector::square_dgemm(int n, const double *A, const double *B, double *C) {
    /* Load data into vectors */
    int adjustN, vN; //Will hold number of elements in a row of the vector version of the output matrix
    Vec4d *vA, *vB, *vC;
    /* vA: n rows x adjustN columns
     * vB: adjustN rows x vN columns
     * vC: n rows x vN columns */

    load_vectors(n, adjustN, vN, A, B, vA, vB, vC);

    /* Perform blocked dgemm using vectors */
    vector_dgemm(n, adjustN, vN, vA, vB, vC);

    /* Extract data from vectors to final matrix */
    store_vectors(n, adjustN, vN, C, vC);

    /* Delete the vector arrays after no longer used */
    delete[] vA;
    delete[] vB;
    delete[] vC;
}

void DgemmVector::load_vectors(int n, int &adjustN, int &vN, const double *A, const double *B, Vec4d *&vA, Vec4d *&vB,
                               Vec4d *&vC) {
    /* Pad the scalar array with zeros if its size is not a multiple of the vector size */
    double *padA, *padB;
    /* padA: n rows x adjustN columns
     * padB: adjustN rows x adjustN columns
     * padC: n rows x adjustN columns */

    /* Pad the matrices with zeros to evenly fill vectors */
    pad_scalar_mats(n, adjustN, A, B, padA, padB);

    /* Load the vector version of A by expanding the cyclic permutation of each chunk of VEC_SIZE adjacent elements of A */
    load_vA(n, adjustN, padA, vA);

    /* Load the vector version of B by looking up elements in diagonal traversal wrapping around VEC_SIZE square blocks of B */
    load_vB(adjustN, vN, padB, vB);

    /* Initialize vC */
    load_vC(n, adjustN, vN, vC);
}

void DgemmVector::store_vectors(int n, int adjustN, int vN, double *C, const Vec4d *vC) {
    double padC[n * adjustN];

    /* Store padded result */
    for(int i=0; i<n; i++)
        for(int j=0; j<vN; j++)
            vC[i + j * n].store(padC + i * VEC_SIZE + j * n);

    /* Unpad result while copying to final result */
    for(int i=0; i<adjustN; i++)
        for(int j=0; j<n; j++)
            C[i + j * n] = padC[i + j * n];
}

void DgemmVector::pad_scalar_mats(int n, int &adjustN, const double *A, const double *B, double *&padA, double *&padB) {
    adjustN = n + n % VEC_SIZE;
    padA = new double[n * adjustN];
    padB = new double[adjustN * adjustN];

    /* Copy data and add padding rows on existing columns */
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            padA[i + j * n] = A[i + j * n];
            padB[i + j * n] = B[i + j * n];
        }

        for(int j=n; j < adjustN; j++)
            padA[i + j * n] = padB[i + j * n] = 0.0;
    }

    /* Add padding columns to B */
    for(int i=n; i < adjustN; i++)
        for(int j=0; j < adjustN; j++)
            padB[i + j * adjustN] = 0.0;
}

void DgemmVector::load_vA(int n, int adjustN, double *padA, Vec4d *&vA) {
    vA = new Vec4d[n * adjustN];
    for(int i=0; i<n; i++)
        for(int j=0; j < adjustN; j += VEC_SIZE)
        {
            Vec4d nextV;
            nextV.load(padA + i + j * n);
            for(int k=0; k<VEC_SIZE; k++)
            {
                vA[i + j * n + k] = nextV;
                nextV = permute4<1,2,3,0>(nextV); //Cycle around one position
            }
        }
}

void DgemmVector::load_vB(int adjustN, int &vN, double *padB, Vec4d *&vB) {
    vN = adjustN / VEC_SIZE; //Initialize vN here
    vB = new Vec4d[adjustN * vN];

    Vec4q lookup_block_row = {0, 1, 2, 3};
    Vec4q lookup_block_col = {0, 1, 2, 3};
    for(int i=0; i < adjustN; i++) //Row (incremented by block)
    {
        for(int j=0; j < vN; j+=VEC_SIZE) //Vector
        {
            for(int k=0; k < VEC_SIZE; k++) //Is the below code correct??
                vB[i + k + j * adjustN] = lookup<VEC_SIZE>((i + lookup_block_row) + (j + lookup_block_col) * adjustN, padB);
            lookup_block_col = permute4<1, 2, 3, 0>(lookup_block_col);
        }

    }


}

void DgemmVector::load_vC(int n, int adjustN, int vN, Vec4d *&vC) {
    vC = new Vec4d[n * vN];
    for(int i=0; i<n; i++)
        for(int j=0; j<vN; j++)
            vC[j + i * vN] = 0.0;
}

