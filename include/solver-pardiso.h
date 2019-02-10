
#ifndef SLSPACK_SOLVER_PARDISO_H
#define SLSPACK_SOLVER_PARDISO_H

#include "data-types.h"
#include "utils.h"

#if USE_PARDISO
int slspack_solver_pardiso(SLSPACK_SOLVER *solver);
#endif

#endif
