
#ifndef SLSPACK_SOLVER_AMG_H
#define SLSPACK_SOLVER_AMG_H

#include "data-types.h"
#include "matrix-utils.h"

#if USE_FASP

#define WITH_UMFPACK 0
#define WITH_PARDISO 0
#define WITH_MUMPS   0

#ifndef __STDC_VERSION__
#define __STDC_VERSION__         199901L
#endif

#include <fasp.h>
#include <fasp_functs.h>

void slspack_solver_amg_create(SLSPACK_SOLVER *s);
void slspack_solver_amg_destroy(SLSPACK_SOLVER *s);

int slspack_solver_amg(SLSPACK_SOLVER *solver);
void slspack_solver_amg_set_pars(SLSPACK_SOLVER *s, AMG_param par);

#endif
#endif
