//
// Created by Orion on 2/22/2022.
//

#include "dgemm_vector.h"
#include "dgemm_utils.h"

void DgemmVector::square_dgemm(int n, const double *A, const double *B, double *C) {
    /* Load data into vectors */
    int adjustN, vN; //Will hold number of elements in a row of the vector version of the output matrix
    Vec4d *vA, *vB, *vC;
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
            vC[j + i * vN].store(padC + j * VEC_SIZE + i * n);

    /* Unpad result while copying to final result */
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            C[j + i * n] = padC[j + i * adjustN];
}

void DgemmVector::pad_scalar_mats(int n, int &adjustN, const double *A, const double *B, double *&padA, double *&padB) {
    adjustN = n + n % VEC_SIZE;
    padA = new double[n * adjustN];
    padB = new double[adjustN * adjustN];

    /* Copy data and add padding columns on existing rows */
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            padA[j + i * adjustN] = A[j + i * n];
            padB[j + i * adjustN] = B[j + i * n];
        }

        for(int j=n; j < adjustN; j++)
            padA[j + i * adjustN] = padB[j + i * adjustN] = 0.0;
    }

    /* Add padding rows */
    for(int i=n; i < adjustN; i++)
        for(int j=0; j < adjustN; j++)
            padB[j + i * adjustN] = 0.0;
}

void DgemmVector::load_vA(int n, int adjustN, double *padA, Vec4d *&vA) {
    vA = new Vec4d[n * adjustN];
    for(int i=0; i<n; i++)
        for(int j=0; j < adjustN; j += VEC_SIZE)
        {
            Vec4d nextV;
            nextV.load(padA + j + i * adjustN);
            for(int k=0; k<VEC_SIZE; k++)
            {
                vA[k + j + i * adjustN] = nextV;
                nextV = permute4<1,2,3,0>(nextV); //Cycle around one position
            }
        }
}

void DgemmVector::load_vB(int adjustN, int &vN, double *padB, Vec4d *&vB) {
    vN = adjustN / VEC_SIZE; //Initialize vN here
    vB = new Vec4d[adjustN * vN];
    for(int i=0; i < adjustN; i += VEC_SIZE)
        for(int j=0; j < vN; j++)
        {
            vB[j + (i) * vN] =     gather4d<0x0, 0x5, 0xa, 0xf>(padB + j * 4 + (i) * adjustN);
            vB[j + (i + 1) * vN] = gather4d<0x4, 0x9, 0xe, 0x3>(padB + j * 4 + (i + 1) * adjustN);
            vB[j + (i + 2) * vN] = gather4d<0x8, 0xd, 0x2, 0x7>(padB + j * 4 + (i + 2) * adjustN);
            vB[j + (i + 3) * vN] = gather4d<0xc, 0x1, 0x6, 0xb>(padB + j * 4 + (i + 3) * adjustN);

            //TODO need to modify if VEC_SIZE is different
        }
}

void DgemmVector::load_vC(int n, int adjustN, int vN, Vec4d *&vC) {
    vC = new Vec4d[n * vN];
    for(int i=0; i<n; i++)
        for(int j=0; j<vN; j++)
            vC[j + i * vN] = 0.0;
}

