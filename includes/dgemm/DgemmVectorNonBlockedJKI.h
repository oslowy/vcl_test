//
// Created by Orion on 2/28/2022.
//

#ifndef VCL_TEST_DGEMMVECTORNONBLOCKEDJKI_H
#define VCL_TEST_DGEMMVECTORNONBLOCKEDJKI_H


#include "DgemmVector.h"

class DgemmVectorNonBlockedJKI: public DgemmVector {
public:
    const char* dgemm_desc() override;
private:
    void vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) override;
};


#endif //VCL_TEST_DGEMMVECTORNONBLOCKEDJKI_H
