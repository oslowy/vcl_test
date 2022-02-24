//
// Created by Serendipity_2 on 2/24/2022.
//

#ifndef VCL_TEST_DGEMMVECTORBLOCKED_H
#define VCL_TEST_DGEMMVECTORBLOCKED_H


#include "dgemm_vector.h"

class DgemmVectorBlocked: public DgemmVector {
protected:
    virtual void do_block(int lda, int M, int N, int K, const Vec4d* A, const Vec4d* B, Vec4d* C) = 0;
};


#endif //VCL_TEST_DGEMMVECTORBLOCKED_H
