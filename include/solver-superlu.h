
#ifndef SLSPACK_SOLVER_SUPERLU_H
#define SLSPACK_SOLVER_SUPERLU_H

#include "data-types.h"
#include "utils.h"

#if USE_SUPERLU

void slspack_solver_superlu_create(SLSPACK_SOLVER *s);
void slspack_solver_superlu_destroy(SLSPACK_SOLVER *s);
int slspack_solver_superlu(SLSPACK_SOLVER *solver);

#endif

#endif
