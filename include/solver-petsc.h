
#ifndef SLSPACK_SOLVER_PETSC_H
#define SLSPACK_SOLVER_PETSC_H

#include "data-types.h"
#include "utils.h"

#if USE_PETSC

typedef void (*solver_petsc_setting)(void *ksp, void *pc);

void slspack_solver_petsc_create(SLSPACK_SOLVER *s);
void slspack_solver_petsc_destroy(SLSPACK_SOLVER *s);

int slspack_solver_petsc(SLSPACK_SOLVER *solver);

void slspack_solver_petsc_setting(solver_petsc_setting func);

#endif

#endif
