
#ifndef SLSPACK_MATRIX_UTILS_H
#define SLSPACK_MATRIX_UTILS_H

#include "data-types.h"
#include "utils.h"

/* free functions */
void slspack_mat_destroy(slspack_mat_csr *A);

/* create csr matrix */
slspack_mat_csr slspack_mat_create(int nrows, int ncols, int *Ap, int *Aj, double *Ax);

/* transpose */
slspack_mat_csr slspack_mat_transpose(const slspack_mat_csr A);

#endif
