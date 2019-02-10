
#ifndef SLSPACK_SOLVER_LIS_H
#define SLSPACK_SOLVER_LIS_H

#include "data-types.h"
#include "utils.h"

#if USE_LIS

void slspack_solver_lis_create(SLSPACK_SOLVER *solver);
void slspack_solver_lis_destroy(SLSPACK_SOLVER *solver);

int slspack_solver_lis(SLSPACK_SOLVER *solver);

void slspack_solver_lis_set_pars(SLSPACK_SOLVER *solver, unsigned int solver_id, unsigned int pc_id);
void slspack_solver_lis_set_option(SLSPACK_SOLVER *solver, char *o);

#endif

#endif
