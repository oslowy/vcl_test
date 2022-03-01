//
// Created by Orion on 2/28/2022.
//

#ifndef VCL_TEST_DGEMMVECTORNONBLOCKEDKIJ_H
#define VCL_TEST_DGEMMVECTORNONBLOCKEDKIJ_H


#include "dgemm_vector.h"

class DgemmVectorNonBlockedKIJ: public DgemmVector {
public:
    const char* dgemm_desc() override;
private:
    void vector_dgemm(int n, int adjustN, int vN, const Vec4d *vA, const Vec4d *vB, Vec4d *vC) override;
};


#endif //VCL_TEST_DGEMMVECTORNONBLOCKEDKIJ_H
