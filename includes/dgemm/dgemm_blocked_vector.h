//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMM_BLOCKED_VECTOR_H
#define VCL_TEST_DGEMM_BLOCKED_VECTOR_H

#include "dgemm_blocked.h"
#include "vectorclass.h"

class DgemmBlockedVector : public DgemmBlocked {
public:
    static const int BLOCK_SIZE = 32;
    static const int VEC_SIZE = 4;
    const char* dgemm_desc() override;

    void square_dgemm (int n, const double* A, const double* B, double* C) override;

private:
    static void load_vectors (int n, const double* A, const double* B, const double* C,
                              int vec_mat_flat, Vec4d* vA, Vec4d* vB, Vec4d* vC);
};

#endif //VCL_TEST_DGEMM_BLOCKED_VECTOR_H
