
#ifndef SLSPACK_SOLVER_MI20_H
#define SLSPACK_SOLVER_MI20_H

#include "data-types.h"
#include "matrix-utils.h"

#if USE_HSL_MI20
int slspack_solver_mi20amg(SLSPACK_SOLVER *solver);
#endif

#endif
