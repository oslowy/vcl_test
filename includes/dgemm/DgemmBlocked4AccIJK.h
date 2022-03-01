//
// Created by Orion on 2/28/2022.
//

#ifndef VCL_TEST_DGEMMBLOCKED4ACCIJK_H
#define VCL_TEST_DGEMMBLOCKED4ACCIJK_H


#include "DgemmBlocked.h"

class DgemmBlocked4AccIJK: public DgemmBlocked {
public:
    const char* dgemm_desc() override;
protected:
    void do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C) override;
};


#endif //VCL_TEST_DGEMMBLOCKED4ACCIJK_H
