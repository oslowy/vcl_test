//
// Created by Serendipity_2 on 2/24/2022.
//

#ifndef VCL_TEST_DGEMMVECTORBLOCKED_H
#define VCL_TEST_DGEMMVECTORBLOCKED_H


#include "DgemmVector.h"

class DgemmVectorBlocked: public DgemmVector {
public:
    static const int BLOCK_SIZE = 8;

protected:
    void vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) override;
    virtual void do_block(int vM, int n, int M, int N, int K, const Vec4d *vA, const double *B, Vec4d *vC) = 0;
};


#endif //VCL_TEST_DGEMMVECTORBLOCKED_H
