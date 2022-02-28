//
// Created by Orion on 2/8/2022.
//

#ifndef VCL_TEST_DGEMM_VECTOR_H
#define VCL_TEST_DGEMM_VECTOR_H

#include "dgemm.h"
#include "vectorclass.h"

class DgemmVector : public Dgemm {
public:
    static const int VEC_SIZE = 4;

    void square_dgemm(int n, const double *A, const double *B, double *C) override;

protected:
    static void load_vectors(int n, int &vN, const double *A, const double *B, Vec4d *&vA, Vec4d *&vB, Vec4d *&vC);
    virtual void vector_dgemm(int n, const Vec4d* vA, const Vec4d* vB, Vec4d* vC) = 0;
    static void store_vectors(int n, double *C, const Vec4d *vC);
};

#endif //VCL_TEST_DGEMM_VECTOR_H
