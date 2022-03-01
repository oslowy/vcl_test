//
// Created by Orion on 2/28/2022.
//

#ifndef VCL_TEST_DGEMMNONBLOCKED4ACCIJK_H
#define VCL_TEST_DGEMMNONBLOCKED4ACCIJK_H


#include "DgemmNonBlocked.h"

class DgemmNonBlocked4AccIJK: public DgemmNonBlocked {
public:
    void square_dgemm (int n, const double* A, const double* B, double* C) override;
    const char* dgemm_desc() override;
};


#endif //VCL_TEST_DGEMMNONBLOCKED4ACCIJK_H
