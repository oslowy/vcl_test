//
// Created by Orion on 2/22/2022.
//

#ifndef VCL_TEST_DGEMMBLOCKED_H
#define VCL_TEST_DGEMMBLOCKED_H

#include "dgemm.h"

class DgemmBlocked : public Dgemm {
public:
    const int BLOCK_SIZE = 32;

protected:
    void square_dgemm (int n, const double* A, const double* B, double* C) override;
    virtual void do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C) = 0;
};

#endif //VCL_TEST_DGEMMBLOCKED_H
