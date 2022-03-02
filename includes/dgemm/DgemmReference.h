//
// Created by Orion on 3/2/2022.
//

#ifndef VCL_TEST_DGEMMREFERENCE_H
#define VCL_TEST_DGEMMREFERENCE_H


#include "Dgemm.h"

class DgemmReference: public Dgemm {
public:
    void square_dgemm (int n, const double* A, const double* B, double* C) override;
    const char* dgemm_desc() override;
};


#endif //VCL_TEST_DGEMMREFERENCE_H
