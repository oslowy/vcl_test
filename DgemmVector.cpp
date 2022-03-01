//
// Created by Orion on 2/22/2022.
//

#include "DgemmVector.h"
#include "dgemm_utils.h"

void DgemmVector::square_dgemm(int n, const double *A, const double *B, double *C) {
    /* Load data into vectors */
    int vN = n + (VEC_SIZE - (n % VEC_SIZE)) % VEC_SIZE; //Round up to a multiple of VEC_SIZE
    int vM = vN / VEC_SIZE; //Will hold number of rows of the vector version of the output matrix
    Vec4d vA[vM * vN], vB[vN * n], vC[vM * n];
    /* vA: vM rows x vN columns
     * vB: vN rows x vK columns
     * vC: vM rows x vK columns */

    load_vectors(vM, vN, n, A, B, vA, vB, vC);

    /* Perform blocked dgemm using vectors */
    vector_dgemm(vM, vN, n, vA, vB, vC);

    /* Extract data from vectors to final matrix */
    store_vectors(vM, vN, n, C, vC);
}

void
DgemmVector::load_vectors(int vM, int vN, int vK, const double *A, const double *B, Vec4d *vA, Vec4d *vB, Vec4d *vC) {
    /* Pad the scalar array with zeros if its size is not a multiple of the vector size */
    double padA[vK * vN], padB[vN * vN];
    /* padA: vK rows x vN columns
     * padB: vN rows x vN columns
     * padC: vK rows x vN columns */

    /* Pad the matrices with zeros to evenly fill vectors */
    pad_scalar_mats(vN, vK, A, B, padA, padB);

    /* Load the vector version of A by expanding the cyclic permutation of each chunk of VEC_SIZE adjacent elements of A */
    load_vA(vM, vN, padA, vA);

    /* Load the vector version of B by looking up elements in diagonal traversal wrapping around VEC_SIZE square blocks of B */
    load_vB(vN, vK, padB, vB);

    /* Initialize vC */
    load_vC(vM, vK, vC);
}

void DgemmVector::store_vectors(int vM, int vN, int vK, double *C, const Vec4d *vC) {
    double padC[vN * vK];

    /* Store padded result */
    for(int i=0; i<vM; i++)
        for(int j=0; j<vK; j++)
            vC[i + j * vM].store(padC + i * VEC_SIZE + j * vN);

    /* Unpad result while copying to final result */
    for(int i=0; i<vK; i++)
        for(int j=0; j<vK; j++)
            C[i + j * vK] = padC[i + j * vN];
}

void DgemmVector::pad_scalar_mats(int vN, int vK, const double *A, const double *B, double *padA, double *padB) {
    /* Copy existing data */
    for(int i=0; i<vK; i++)
        for(int j=0; j<vK; j++)
        {
            padA[i + j * vN] = A[i + j * vK];
            padB[i + j * vN] = B[i + j * vK];
        }

    /* Add padding rows to A and B and padding columns to A */
    for(int i=vK; i < vN; i++)
    {
        for(int j=0; j < vK; j++)
            padA[i + j * vN] = padB[i + j * vN] = 0.0;

        for(int j=0; j < vN; j++)
            padA[i + j * vN] = 0.0;
    }
}

void DgemmVector::load_vA(int vM, int vN, double *padA, Vec4d *vA) {
    Vec4q lookup_block_row = {0, 1, 2, 3};
    Vec4q lookup_block_col = {0, 1, 2, 3};

    for(int i=0; i < vM; i++) //Row of vA
        for(int j=0; j < vN; j += VEC_SIZE) //Column of vA (incremented by block of VEC_SIZE)
        {
            /* Scramble the vectors so that the elements will multiply with the permuted
             * copies of vB in the correct order */
            for(int k=0; k < VEC_SIZE; k++)
            {
                Vec4q lookup_i = i * VEC_SIZE + lookup_block_row + (j + lookup_block_col) * vN;
                Vec4d next_v = lookup<VEC_SIZE>(lookup_i, padA);
                vA[i + (j + k) * vM] = next_v;
                lookup_block_col = permute4<1, 2, 3, 0>(lookup_block_col);
            }
        }
}

void DgemmVector::load_vB(int vN, int vK, double *padB, Vec4d *vB) {

    for(int i=0; i < vN; i += VEC_SIZE) //Row (incremented by block of VEC_SIZE)
    {
        for(int j=0; j < vK; j++) //Column
        {
            /* Copy each chunk of four elements around four permutations of 1 cyclic shift so that each element will
             * participate in a vector multiply with each corresponding element in B*/
            Vec4d nextV;
            nextV.load(padB + i + j * vN);
            for(int k=0; k<VEC_SIZE; k++)
            {
                vB[i + k + j * vN] = nextV;
                nextV = permute4<1, 2, 3, 0>(nextV); //Cycle around one position
            }
        }
    }
}

void DgemmVector::load_vC(int vM, int vK, Vec4d *vC) {
    for(int i=0; i<vM; i++)
        for(int j=0; j<vK; j++)
            vC[j + i * vM] = 0.0;
}

