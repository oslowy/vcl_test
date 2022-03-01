//
// Created by Orion on 2/28/2022.
//

#ifndef VCL_TEST_DGEMMVECTORBLOCKEDJKI_H
#define VCL_TEST_DGEMMVECTORBLOCKEDJKI_H

#include "DgemmVectorBlocked.h"

class DgemmVectorBlockedJKI : public DgemmVectorBlocked {
public:
    const char* dgemm_desc() override;
protected:
    void do_block(int vM, int vN, int vK, int M, int N, int K, const Vec4d *vA, const Vec4d *vB, Vec4d *vC) override;
};


#endif //VCL_TEST_DGEMMVECTORBLOCKEDJKI_H
