//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMMVECTOR_H
#define VCL_TEST_DGEMMVECTOR_H

#include "Dgemm.h"
#include "vectorclass.h"

#ifndef VEC_SIZE
#define VEC_SIZE 4
#endif

void pad_scalar_mat(int padN, int n, double *padA, const double *A);
void load_vA(int vM, int padN, int n, Vec4d *vA, double *padA);
void load_vC(int vM, int n, Vec4d *vC);

void load_vectors(int vM, int n, const double *A, Vec4d *vA, Vec4d *vC);
void vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC);
void store_vectors(int vM, int n, double *C, const Vec4d *vC);

#endif //VCL_TEST_DGEMMVECTOR_H
