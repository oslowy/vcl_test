//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMM_BLOCKED_VECTOR_H
#define VCL_TEST_DGEMM_BLOCKED_VECTOR_H

#include "dgemm_blocked.h"

class DgemmBlockedVector : public DgemmBlocked {
public:
    static const int BLOCK_SIZE = 32;
    static const int VEC_SIZE = 4;
    const char* dgemm_desc() override;

    void do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C) override;
};

#endif //VCL_TEST_DGEMM_BLOCKED_VECTOR_H
