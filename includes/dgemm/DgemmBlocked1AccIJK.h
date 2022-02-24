//
// Created by Serendipity_2 on 2/24/2022.
//

#ifndef VCL_TEST_DGEMMBLOCKED1ACCIJK_H
#define VCL_TEST_DGEMMBLOCKED1ACCIJK_H


#include "dgemm_blocked.h"

class DgemmBlocked1AccIJK: public DgemmBlocked {
    void do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C) override;
    void square_dgemm (int n, const double* A, const double* B, double* C) override;
    const char* dgemm_desc() override;
};


#endif //VCL_TEST_DGEMMBLOCKED1ACCIJK_H
