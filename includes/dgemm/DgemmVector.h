//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMMVECTOR_H
#define VCL_TEST_DGEMMVECTOR_H

#include "Dgemm.h"
#include "vectorclass.h"

class DgemmVector : public Dgemm {
public:
    static const int VEC_SIZE = 4;

    void square_dgemm(int n, const double *A, const double *B, double *C) override;

protected:
    static void pad_scalar_mat(int padN, int n, double *padA, const double *A);
    static void load_vA(int vM, int n, Vec4d *vA, double *padA);
    static void load_vC(int vM, int n, Vec4d *vC);

    static void load_vectors(int vM, int n, const double *A, Vec4d *vA, Vec4d *vC);
    virtual void vector_dgemm(int vM, int n, const Vec4d *vA, const double *B, Vec4d *vC) = 0;
    static void store_vectors(int vM, int n, double *C, const Vec4d *vC);
};

#endif //VCL_TEST_DGEMMVECTOR_H
