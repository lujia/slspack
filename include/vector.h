
#ifndef SLSPACK_VECTOR_H
#define SLSPACK_VECTOR_H

#include "matrix-utils.h"

/* create and destroy vector */
slspack_vec slspack_vec_create(int n);
void slspack_vec_destroy(slspack_vec *v);

/* add value */
void slspack_vec_add_value_by_index(slspack_vec x, int i, double val);

/* set value */
void slspack_vec_set_value(slspack_vec x, double val);
void slspack_vec_set_value_by_array(slspack_vec x, double *val);
void slspack_vec_set_value_by_index(slspack_vec x, int i, double val);

/* get value from vector */
void slspack_vec_get_value(double *val, slspack_vec x);
double slspack_vec_get_value_by_index(slspack_vec x, int i);

/* copy */
void slspack_vec_copy(slspack_vec des, const slspack_vec src);

double slspack_vec_norm(slspack_vec x);

#endif
