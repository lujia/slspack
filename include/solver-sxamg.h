
#ifndef SLSPACK_SOLVER_SXAMG_H
#define SLSPACK_SOLVER_SXAMG_H

#include "data-types.h"
#include "matrix-utils.h"

#if USE_SXAMG

#include "sxamg.h"

void slspack_solver_sxamg_create(SLSPACK_SOLVER *s);
void slspack_solver_sxamg_destroy(SLSPACK_SOLVER *s);

int slspack_solver_sxamg(SLSPACK_SOLVER *solver);

void slspack_solver_sxamg_set_pars(SLSPACK_SOLVER *solver, SX_AMG_PARS *pars);

#endif

#endif
