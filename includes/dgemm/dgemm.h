//
// Created by Orion on 2/22/2022.
//

#ifndef VCL_TEST_DGEMM_H
#define VCL_TEST_DGEMM_H


class Dgemm {
public:
    virtual void square_dgemm (int n, const double* A, const double* B, double* C) = 0;
    const char* dgemm_desc = "";
};


#endif //VCL_TEST_DGEMM_H
