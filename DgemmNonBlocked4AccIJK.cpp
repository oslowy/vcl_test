//
// Created by Orion on 2/28/2022.
//

#include "DgemmNonBlocked4AccIJK.h"

void DgemmNonBlocked4AccIJK::square_dgemm(int n, const double *A, const double *B, double *C) {

}

const char *DgemmNonBlocked4AccIJK::dgemm_desc() {
    return "Non-blocked DGEMM with i-j-k loop order and four accumulators";
}
