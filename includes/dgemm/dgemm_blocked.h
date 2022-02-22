//
// Created by Orion on 2/22/2022.
//

#ifndef VCL_TEST_DGEMM_BLOCKED_H
#define VCL_TEST_DGEMM_BLOCKED_H

#include "dgemm.h"

class DgemmBlocked : public Dgemm {
public:
    static const int BLOCK_SIZE = 32;

    virtual void do_block (int lda, int M, int N, int K, double* A, double* B, double* C);
    void square_dgemm (int lda, double* A, double* B, double* C);
};

#endif //VCL_TEST_DGEMM_BLOCKED_H
