
#ifndef SLSPACK_SOLVER_FASP_H
#define SLSPACK_SOLVER_FASP_H

#include "data-types.h"
#include "utils.h"

#if USE_FASP

#define WITH_UMFPACK             0
#define WITH_PARDISO             0
#define WITH_MUMPS               0

#include <fasp.h>
#include <fasp_functs.h>

void slspack_solver_fasp_create(SLSPACK_SOLVER *s);
void slspack_solver_fasp_destroy(SLSPACK_SOLVER *s);

void slspack_solver_fasp_set_pars(SLSPACK_SOLVER *solver, input_param *inparam, ITS_param *itsparam,
        AMG_param *amgparam, ILU_param *iluparam, SWZ_param *schparam);

int slspack_solver_fasp(SLSPACK_SOLVER *solver);

#endif

#endif
