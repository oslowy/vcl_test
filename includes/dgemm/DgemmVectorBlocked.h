//
// Created by Serendipity_2 on 2/24/2022.
//

#ifndef VCL_TEST_DGEMMVECTORBLOCKED_H
#define VCL_TEST_DGEMMVECTORBLOCKED_H


#include "DgemmVector.h"

#ifndef VEC_BLOCK_SIZE
#define VEC_BLOCK_SIZE 8
#endif

void do_block(int vM, int n, int M, int N, int K, const Vec4d *vA, const double *B, Vec4d *vC);

#endif //VCL_TEST_DGEMMVECTORBLOCKED_H
