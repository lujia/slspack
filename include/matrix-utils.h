
#ifndef SLSPACK_MATRIX_UTILS_H
#define SLSPACK_MATRIX_UTILS_H

#include "data-types.h"
#include "utils.h"

/* free functions */
void slspack_mat_destroy(SLSPACK_MAT *A);

/* create csr matrix */
SLSPACK_MAT slspack_mat_create(int nrows, int ncols, int *Ap, int *Aj, double *Ax);

/* transpose */
SLSPACK_MAT slspack_mat_transpose(const SLSPACK_MAT A);

#endif
