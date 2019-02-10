
#ifndef SLSPACK_SOLVER_LASPACK_H
#define SLSPACK_SOLVER_LASPACK_H

#include "data-types.h"
#include "utils.h"

#if USE_LASPACK
void slspack_solver_laspack_create(SLSPACK_SOLVER *s);
void slspack_solver_laspack_destroy(SLSPACK_SOLVER *s);

int slspack_solver_laspack(SLSPACK_SOLVER *solver);

void slspack_solver_laspack_set_pars(SLSPACK_SOLVER *solver, int solver_type, int pc_type);
#endif

#endif
