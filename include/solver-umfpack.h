
#ifndef SLSPACK_SOLVER_UMFPACK_H
#define SLSPACK_SOLVER_UMFPACK_H

#include "data-types.h"
#include "utils.h"

#if USE_SSPARSE

void slspack_solver_umfpack_create(SLSPACK_SOLVER *s);
void slspack_solver_umfpack_destroy(SLSPACK_SOLVER *s);

int slspack_solver_umfpack(SLSPACK_SOLVER *solver);

#endif

#endif
