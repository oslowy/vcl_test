//
// Created by Orion on 2/28/2022.
//

#ifndef VCL_TEST_DGEMMNONBLOCKED4ACCJKI_H
#define VCL_TEST_DGEMMNONBLOCKED4ACCJKI_H


#include "DgemmNonBlocked.h"

class DgemmNonBlocked4AccJKI: public DgemmNonBlocked {
public:
    void square_dgemm (int n, const double* A, const double* B, double* C) override;
    const char* dgemm_desc() override;
};


#endif //VCL_TEST_DGEMMNONBLOCKED4ACCJKI_H
