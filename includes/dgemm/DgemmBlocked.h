//
// Created by Orion on 2/22/2022.
//

#ifndef VCL_TEST_DGEMMBLOCKED_H
#define VCL_TEST_DGEMMBLOCKED_H

#include "Dgemm.h"

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 32
#endif

void do_block (int lda, int M, int N, int K, const double* A, const double* B, double* C);

#endif //VCL_TEST_DGEMMBLOCKED_H
