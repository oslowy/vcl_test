//
// Created by Orion on 2/28/2022.
//

#include "DgemmNonBlocked4AccJKI.h"

void DgemmNonBlocked4AccJKI::square_dgemm(int n, const double *A, const double *B, double *C) {

}

const char *DgemmNonBlocked4AccJKI::dgemm_desc() {
    return "Non-blocked DGEMM with j-k-i loop order and four accumulators";
}
