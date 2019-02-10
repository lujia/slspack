
#ifndef SLSPACK_SOLVER_MUMPS_H
#define SLSPACK_SOLVER_MUMPS_H

#include "data-types.h"
#include "utils.h"

#if USE_MUMPS
void slspack_solver_mumps_create(SLSPACK_SOLVER *s);
void slspack_solver_mumps_destroy(SLSPACK_SOLVER *s);

int slspack_solver_mumps(SLSPACK_SOLVER *solver);
#endif

#endif
