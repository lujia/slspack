
#ifndef SLSPACK_SLSPACK_H
#define SLSPACK_SLSPACK_H

#include "vector.h"
#include "solver-laspack.h"
#include "solver-mumps.h"
#include "solver-petsc.h"
#include "solver-lis.h"
#include "solver-umfpack.h"
#include "solver-klu.h"
#include "solver-fasp.h"
#include "solver-superlu.h"
#include "solver-pardiso.h"
#include "solver-mi20.h"
#include "solver-amg.h"
#include "solver-sxamg.h"

/* create solver */
void slspack_solver_create(SLSPACK_SOLVER *s, SLSPACK_SOLVER_TYPE s_type);

/* assemble solver */
void slspack_solver_assemble(SLSPACK_SOLVER *s, slspack_mat_csr A, slspack_vec x, slspack_vec b);

/* destroy solver */
void slspack_solver_destroy(SLSPACK_SOLVER *s);

/* solve */
int slspack_solver_solve(SLSPACK_SOLVER *solver);

/* set relative tolerence */
void slspack_solver_set_rtol(SLSPACK_SOLVER *s, double tol);

/* set absolute tolerence */
void slspack_solver_set_atol(SLSPACK_SOLVER *s, double tol);

/* set relative b norm tolerence */
void slspack_solver_set_rbtol(SLSPACK_SOLVER *s, double tol);

/* set maximal number of iteration */
void slspack_solver_set_maxit(SLSPACK_SOLVER *s, int maxit);

/* set the number of restart */
void slspack_solver_set_restart(SLSPACK_SOLVER *s, int m);

/* reset verbosity of solver */
void slspack_solver_set_verbosity(SLSPACK_SOLVER *s, int v);

double slspack_solver_get_residual(SLSPACK_SOLVER s);
int slspack_solver_get_nits(SLSPACK_SOLVER s);

#endif
