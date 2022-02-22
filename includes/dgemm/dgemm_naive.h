//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMM_NAIVE_H
#define VCL_TEST_DGEMM_NAIVE_H

#include "dgemm_utils.h"
#include "dgemm.h"

class DgemmNaive : public Dgemm {
public:
    void square_dgemm (int n, const double* A, const double* B, double* C) override;
    const char* dgemm_desc {"Naive DGEMM"};
};

#endif //VCL_TEST_DGEMM_NAIVE_H
