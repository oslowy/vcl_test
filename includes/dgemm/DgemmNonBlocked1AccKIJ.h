//
// Created by Serendipity_2 on 2/24/2022.
//

#ifndef VCL_TEST_DGEMMNONBLOCKED1ACCKIJ_H
#define VCL_TEST_DGEMMNONBLOCKED1ACCKIJ_H

#include "dgemm_utils.h"
#include "dgemm.h"
#include "dgemm_naive.h"

class DgemmNonBlocked1AccKIJ: public DgemmNonBlocked {
public:
    void square_dgemm (int n, const double* A, const double* B, double* C) override;
    const char* dgemm_desc() override;
};


#endif //VCL_TEST_DGEMMNONBLOCKED1ACCKIJ_H