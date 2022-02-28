//
// Created by Orion on 2/22/2022.
//

#include "dgemm_vector.h"
#include "dgemm_utils.h"

void DgemmVector::square_dgemm(int n, const double *A, const double *B, double *C) {
    /* Load data into vectors */
    int vN; //Will hold number of elements in the vector version of the matrix
    Vec4d *vA, *vB, *vC;
    load_vectors(n, vN, A, B, vA, vB, vC);

    /* Perform blocked dgemm using vectors */
    vector_dgemm(n, vA, vB, vC);

    /* Extract data from vectors to final matrix */
    store_vectors(n, C, vC);

    /* Delete the vector arrays after no longer used */
    delete[] vA;
    delete[] vB;
    delete[] vC;
}

void DgemmVector::load_vectors(int n, int &vN, const double *A, const double *B, Vec4d *&vA, Vec4d *&vB, Vec4d *&vC) {
    /* Pad the scalar array with zeros if its size is not a multiple of the vector size */
    int remainder = n % VEC_SIZE;
    vN = n + remainder;
    double padA[n * vN], padB[vN * vN];

    /* Copy data and add padding columns on existing rows */
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            padA[j + i * vN] = A[j + i * n];
            padB[j + i * vN] = B[j + i * n];
        }

        for(int j=n; j < vN; j++)
            padA[j + i * vN] = padB[j + i * vN] = 0.0;
    }

    /* Add padding rows */
    for(int i=n; i < vN; i++)
        for(int j=0; j < vN; j++)
            padA[j + i * vN] = padB[j + i * vN] = 0.0;

    /* Load the vector version of A by expanding the cyclic permutation of each chunk of VEC_SIZE adjacent elements of A */
    vA = new Vec4d[n * vN];
    for(int i=0; i<n; i++)
        for(int j=0; j < vN; j += VEC_SIZE)
        {
            Vec4d nextV;
            nextV.load(padA + j + i * vN);
            for(int k=0; k<VEC_SIZE; k++)
            {
                vA[k + j + i * vN] = nextV;
                nextV = permute4<1,2,3,0>(nextV); //Cycle around one position
            }
        }

    /* Load the vector version of B using a special lookup pattern */
    vB = new Vec4d[vN * vN];


}

void DgemmVector::store_vectors(int n, double *C, const Vec4d *vC) {

}
