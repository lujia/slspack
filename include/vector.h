
#ifndef SLSPACK_VECTOR_H
#define SLSPACK_VECTOR_H

#include "matrix-utils.h"

/* create and destroy vector */
SLSPACK_VEC slspack_vec_create(int n);
void slspack_vec_destroy(SLSPACK_VEC *v);

/* add value */
void slspack_vec_add_value_by_index(SLSPACK_VEC x, int i, double val);

/* set value */
void slspack_vec_set_value(SLSPACK_VEC x, double val);
void slspack_vec_set_value_by_array(SLSPACK_VEC x, double *val);
void slspack_vec_set_value_by_index(SLSPACK_VEC x, int i, double val);

/* get value from vector */
void slspack_vec_get_value(double *val, SLSPACK_VEC x);
double slspack_vec_get_value_by_index(SLSPACK_VEC x, int i);

/* copy */
void slspack_vec_copy(SLSPACK_VEC des, const SLSPACK_VEC src);

double slspack_vec_norm(SLSPACK_VEC x);

#endif
