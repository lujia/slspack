
#ifndef SLSPACK_SOLVER_KLU_H
#define SLSPACK_SOLVER_KLU_H

#include "data-types.h"
#include "utils.h"
#include "matrix-utils.h"

#if USE_SSPARSE
int slspack_solver_klu(SLSPACK_SOLVER *solver);
#endif

#endif
