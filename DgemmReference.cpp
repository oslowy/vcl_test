//
// Created by Orion on 3/2/2022.
//

#include "DgemmReference.h"
#include <cblas.h> // For: cblas_dgemm

void DgemmReference::square_dgemm(int n, const double *A, const double *B, double *C) {
    int M = n;
    int N = n;
    int K = n;
    double ALPHA = 1.;
    double BETA = 1.;
    int LDA = N;
    int LDB = N;
    int LDC = N;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

const char *DgemmReference::dgemm_desc() {
    return "Reference to default DGEMM for comparison";
}
