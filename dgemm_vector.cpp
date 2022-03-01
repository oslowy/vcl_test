//
// Created by Orion on 2/22/2022.
//

#include "dgemm_vector.h"
#include "dgemm_utils.h"

void DgemmVector::square_dgemm(int n, const double *A, const double *B, double *C) {
    /* Load data into vectors */
    int adjustN = n + (VEC_SIZE - (n % VEC_SIZE)) % VEC_SIZE; //Round up to a multiple of VEC_SIZE
    int vN = adjustN / VEC_SIZE; //Will hold number of rows of the vector version of the output matrix
    Vec4d vA[vN * adjustN], vB[adjustN * n], vC[vN * n];
    /* vA: vN rows x adjustN columns
     * vB: adjustN rows x n columns
     * vC: vN rows x n columns */

    load_vectors(n, adjustN, vN, A, B, vA, vB, vC);

    /* Perform blocked dgemm using vectors */
    vector_dgemm(n, adjustN, vN, vA, vB, vC);

    /* Extract data from vectors to final matrix */
    store_vectors(n, adjustN, vN, C, vC);
}

void DgemmVector::load_vectors(int n, int adjustN, int vN, const double *A, const double *B, Vec4d *vA, Vec4d *vB,
                               Vec4d *vC) {
    /* Pad the scalar array with zeros if its size is not a multiple of the vector size */
    double padA[n * adjustN], padB[adjustN * adjustN];
    /* padA: n rows x adjustN columns
     * padB: adjustN rows x adjustN columns
     * padC: n rows x adjustN columns */

    /* Pad the matrices with zeros to evenly fill vectors */
    pad_scalar_mats(n, adjustN, A, B, padA, padB);

    /* Load the vector version of A by expanding the cyclic permutation of each chunk of VEC_SIZE adjacent elements of A */
    load_vA(vN, adjustN, padA, vA);

    /* Load the vector version of B by looking up elements in diagonal traversal wrapping around VEC_SIZE square blocks of B */
    load_vB(adjustN, n, padB, vB);

    /* Initialize vC */
    load_vC(n, vN, vC);
}

void DgemmVector::store_vectors(int n, int adjustN, int vN, double *C, const Vec4d *vC) {
    double padC[adjustN * n];

    /* Store padded result */
    for(int i=0; i<vN; i++)
        for(int j=0; j<n; j++)
            vC[i + j * vN].store(padC + i * VEC_SIZE + j * adjustN);

    /* Unpad result while copying to final result */
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            C[i + j * n] = padC[i + j * adjustN];
}

void DgemmVector::pad_scalar_mats(int n, int &adjustN, const double *A, const double *B, double *padA, double *padB) {
    /* Copy existing data */
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
        {
            padA[i + j * adjustN] = A[i + j * n];
            padB[i + j * adjustN] = B[i + j * n];
        }

    /* Add padding rows to A and B and padding columns to A */
    for(int i=n; i < adjustN; i++)
    {
        for(int j=0; j < n; j++)
            padA[i + j * adjustN] = padB[i + j * adjustN] = 0.0;

        for(int j=0; j < adjustN; j++)
            padA[i + j * adjustN] = 0.0;
    }
}

void DgemmVector::load_vA(int &vN, int adjustN, double *padA, Vec4d *&vA) {
    vA = new Vec4d[vN * adjustN];
    Vec4q lookup_block_row = {0, 1, 2, 3};
    Vec4q lookup_block_col = {0, 1, 2, 3};

    for(int i=0; i < vN; i++) //Row of vA
        for(int j=0; j < adjustN; j += VEC_SIZE) //Column of vA (incremented by block of VEC_SIZE)
        {
            /* Scramble the vectors so that the elements will multiply with the permuted
             * copies of vB in the correct order */
            for(int k=0; k < VEC_SIZE; k++)
            {
                Vec4q lookup_i = i * VEC_SIZE + lookup_block_row + (j + lookup_block_col) * adjustN;
                Vec4d next_v = lookup<VEC_SIZE>(lookup_i, padA);
                vA[i + (j + k) * vN] = next_v;
                lookup_block_col = permute4<1, 2, 3, 0>(lookup_block_col);
            }
        }
}

void DgemmVector::load_vB(int adjustN, int n, double *padB, Vec4d *&vB) {
    vB = new Vec4d[adjustN * n];

    for(int i=0; i < adjustN; i += VEC_SIZE) //Row (incremented by block of VEC_SIZE)
    {
        for(int j=0; j < n; j++) //Column
        {
            /* Copy each chunk of four elements around four permutations of 1 cyclic shift so that each element will
             * participate in a vector multiply with each corresponding element in B*/
            Vec4d nextV;
            nextV.load(padB + i + j * adjustN);
            for(int k=0; k<VEC_SIZE; k++)
            {
                vB[i + k + j * adjustN] = nextV;
                nextV = permute4<1, 2, 3, 0>(nextV); //Cycle around one position
            }
        }
    }
}

void DgemmVector::load_vC(int n, int vN, Vec4d *&vC) {
    vC = new Vec4d[vN * n];
    for(int i=0; i<vN; i++)
        for(int j=0; j<n; j++)
            vC[j + i * vN] = 0.0;
}

