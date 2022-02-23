//
// Created by Orion on 2/22/2022.
//

#include "dgemm_blocked_vector.h"
#include "dgemm_utils.h"

const char *DgemmBlockedVector::dgemm_desc() {
    return "Vectorized Blocked DGEMM";
}

void DgemmBlockedVector::square_dgemm(int n, const double *A, const double *B, double *C) {
    /* Load data into vectors */
    int vec_mat_stop = n / VEC_SIZE; //Length of one row of vector matrix if source matrix dimension is a multiple of vector size
    int vec_mat_dim = (n % VEC_SIZE == 0)? vec_mat_stop: vec_mat_stop + 1;
    int vec_mat_flat = n * vec_mat_dim;

    Vec4d vA[vec_mat_flat], vB[vec_mat_flat], vC[vec_mat_flat];
    load_vectors(n, A, B, C,
                 vec_mat_flat, vA, vB, vC);

    /* Perform blocked dgemm using vectors */

    /* Extract data from vectors to final matrix */

    //Stop reading before padding in padded vectors
}

void DgemmBlockedVector::load_vectors(int n, const double *A, const double *B, const double *C,
                                      int vec_mat_flat, Vec4d *vA, Vec4d *vB, Vec4d *vC) {
    int vec_mat_stop = n / VEC_SIZE;
    int vec_flat_index = 0;
    int scalar_flat_offset = 0;
    int remainder = n % VEC_SIZE;

    for (int i=0; i<n; i++)
    {
        /* Fill main portion of row */
        for (int j=0; j<vec_mat_stop; j++)
        {
            vA[vec_flat_index].load(A + scalar_flat_offset);
            vB[vec_flat_index].load(B + scalar_flat_offset);
            vC[vec_flat_index] = 0.0; //Initialize accumulator

            /* Update index */
            vec_flat_index++;
            scalar_flat_offset += VEC_SIZE;
        }

        /* Fill end of row that partially oversteps original row end */
        vA[vec_flat_index].load_partial(remainder, A + scalar_flat_offset);
        vB[vec_flat_index].load_partial(remainder, B + scalar_flat_offset);
        vC[vec_flat_index] = 0.0;

        /* Update index */
        vec_flat_index++;
        scalar_flat_offset += remainder;
    }

}
