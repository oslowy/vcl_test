//
// Created by Orion on 2/22/2022.
//

#ifndef VCL_TEST_DGEMM_BLOCKED_H
#define VCL_TEST_DGEMM_BLOCKED_H

#include "dgemm.h"
#include "dgemm_utils.h"

class DgemmBlocked : public Dgemm {
public:
    const int BLOCK_SIZE = 32;

    virtual void do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C) = 0;
};

#endif //VCL_TEST_DGEMM_BLOCKED_H
