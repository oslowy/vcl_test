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
    static void pad_scalar_mats(int vN, int vK, const double *A, const double *B, double *padA, double *padB);
    static void load_vA(int vM, int vN, double *padA, Vec4d *vA);
    static void load_vB(int vN, int vK, double *padB, Vec4d *vB);
    static void load_vC(int vM, int vK, Vec4d *vC);

    static void load_vectors(int vM, int vN, int vK, const double *A, const double *B, Vec4d *vA, Vec4d *vB, Vec4d *vC);
    virtual void vector_dgemm(int vM, int vN, int vK, const Vec4d *vA, const Vec4d *vB, Vec4d *vC) = 0;
    static void store_vectors(int vM, int vN, int vK, double *C, const Vec4d *vC);
};

#endif //VCL_TEST_DGEMMVECTOR_H
