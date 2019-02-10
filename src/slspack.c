
#include "slspack.h"

/* constants */
int SLSPACK_RESTART        = 30;
int SLSPACK_MAXIT          = 1000;

double SLSPACK_ATOL        = 1e-7;
double SLSPACK_RTOL        = 1e-7;
double SLSPACK_RBN         = 1e-7;
double SLSPACK_BREAKDOWN   = 1e-40;

void slspack_solver_create(SLSPACK_SOLVER *s, SLSPACK_SOLVER_TYPE s_type)
{
    assert(s != NULL);

    bzero(s, sizeof(SLSPACK_SOLVER));

    s->type = s_type;
    s->residual = 0;

    /* pars */
    s->tol_rel = SLSPACK_RTOL;
    s->tol_abs = SLSPACK_ATOL;
    s->tol_rbn = SLSPACK_RBN;

    s->restart = SLSPACK_RESTART;
    s->maxit = SLSPACK_MAXIT;
    s->nits = 0;

    /* default pars */
#if USE_FASP
    s->fasp_solver = 1;
    s->fasp_pc = 4;
#endif

#if USE_LIS
    s->lis_solver = 0;
    s->lis_pc = 2;
#endif


#if USE_PETSC
    s->petsc_solver = 1;
    s->petsc_pc = 1;
#endif

#if USE_LASPACK
    s->laspack_solver = 0;
    s->laspack_pc = 0;
#endif

    bzero(&s->A, sizeof(s->A));

    s->verb = 2;

#if USE_LASPACK
    if (s_type == SLSPACK_SOLVER_LASPACK) {
        slspack_solver_laspack_create(s);
    }
#endif

#if USE_SSPARSE
    if (s_type == SLSPACK_SOLVER_UMFPACK) {
        slspack_solver_umfpack_create(s);
    }
#endif

#if USE_MUMPS
    if (s_type == SLSPACK_SOLVER_MUMPS) {
        slspack_solver_mumps_create(s);
    }
#endif

#if USE_PETSC
    if (s_type == SLSPACK_SOLVER_PETSC) {
        slspack_solver_petsc_create(s);
    }
#endif

    /* fasp */
#if USE_FASP
    if (s_type == SLSPACK_SOLVER_FASP) {
        slspack_solver_fasp_create(s);
    }

    /* amg */
    if (s_type == SLSPACK_SOLVER_AMG || s_type == SLSPACK_SOLVER_FMG) {
        slspack_solver_amg_create(s);
    }
#endif

#if USE_SUPERLU
    if (s_type == SLSPACK_SOLVER_SUPERLU) {
        slspack_solver_superlu_create(s);
    }
#endif

#if USE_SXAMG
    if (s_type == SLSPACK_SOLVER_SXAMG) {
        slspack_solver_sxamg_create(s);
    }
#endif

#if USE_LIS
    if (s_type == SLSPACK_SOLVER_LIS) {
        slspack_solver_lis_create(s);
    }
#endif

    s->assembled = 0;
}

void slspack_solver_assemble(SLSPACK_SOLVER *s, slspack_mat_csr Ax, slspack_vec x, slspack_vec b)
{
    slspack_mat_csr A;
    double t = 0.;

    assert(s != NULL);

    if (Ax.num_rows <= 0) {
        slspack_error(1, "solver: wrong input matrix, number of rows "
                "should be greater than 0\n");
    }
    else {
        if (Ax.num_rows != Ax.num_cols) {
            slspack_error(1, "solver: wrong input matrix, number of rows != "
                    "number of columns\n");
        }
    }

    if (Ax.num_nnzs < Ax.num_rows) {
        slspack_error(1, "solver: wrong input matrix, singular\n");
    }

    if (s->verb > 1) {
        t = slspack_get_time();
    }

    A.num_rows = Ax.num_rows;
    A.num_cols = Ax.num_cols;
    A.num_nnzs = Ax.num_nnzs;
    A.Ap = slspack_copy_on(Ax.Ap, sizeof(int) * (Ax.num_rows + 1));
    A.Aj = slspack_copy_on(Ax.Aj, sizeof(int) * Ax.num_nnzs);
    A.Ax = slspack_copy_on(Ax.Ax, sizeof(double) * Ax.num_nnzs);

    s->rhs = b;
    s->x = x;

    s->A = A;

    s->assembled = 1;

    if (s->verb > 1) {
        t = slspack_get_time() - t;

        slspack_printf("solver: assemble time: %g\n", t);
    }
}

void slspack_solver_destroy(SLSPACK_SOLVER *s)
{
    assert(s != NULL);

    assert(s->assembled);
    slspack_mat_destroy(&s->A);

    /* fasp */
#if USE_FASP
    if (s->type == SLSPACK_SOLVER_FASP) {
        slspack_solver_fasp_destroy(s);
    }

    /* amg solver */
    if (s->type == SLSPACK_SOLVER_AMG || s->type == SLSPACK_SOLVER_FMG) {
        slspack_solver_amg_destroy(s);
    }
#endif

    /* laspack */
#if USE_LASPACK
    if (s->type == SLSPACK_SOLVER_LASPACK) {
        slspack_solver_laspack_destroy(s);
    }
#endif

    /* mumps */
#if USE_MUMPS
    if (s->type == SLSPACK_SOLVER_MUMPS) {
        slspack_solver_mumps_destroy(s);
    }
#endif

#if USE_PETSC
    if (s->type == SLSPACK_SOLVER_PETSC) {
        slspack_solver_petsc_destroy(s);
    }
#endif

#if USE_SUPERLU
    if (s->type == SLSPACK_SOLVER_SUPERLU) {
        slspack_solver_superlu_destroy(s);
    }
#endif

#if USE_SSPARSE
    if (s->type == SLSPACK_SOLVER_UMFPACK) {
        slspack_solver_umfpack_destroy(s);
    }
#endif

#if USE_SXAMG
    if (s->type == SLSPACK_SOLVER_SXAMG) {
        slspack_solver_sxamg_destroy(s);
    }
#endif

#if USE_LIS
    if (s->type == SLSPACK_SOLVER_LIS) {
        slspack_solver_lis_destroy(s);
    }
#endif

    s->assembled = 0;
}

int slspack_solver_solve(SLSPACK_SOLVER *solver)
{
    SLSPACK_SOLVER_TYPE type = solver->type;

    assert(solver->assembled);

    switch (type) {
#if USE_LASPACK
        case SLSPACK_SOLVER_LASPACK:
            return slspack_solver_laspack(solver);
#endif

#if USE_SSPARSE
        case SLSPACK_SOLVER_UMFPACK:
            return slspack_solver_umfpack(solver);

        case SLSPACK_SOLVER_KLU:
            return slspack_solver_klu(solver);
#endif

#if USE_MUMPS
        case SLSPACK_SOLVER_MUMPS:
            return slspack_solver_mumps(solver);
#endif

#if USE_PETSC
        case SLSPACK_SOLVER_PETSC:
            return slspack_solver_petsc(solver);
#endif

#if USE_LIS
        case SLSPACK_SOLVER_LIS:
            return slspack_solver_lis(solver);
#endif

#if USE_FASP
        case SLSPACK_SOLVER_FASP:
            return slspack_solver_fasp(solver);

        case SLSPACK_SOLVER_AMG:
        case SLSPACK_SOLVER_FMG:
            return slspack_solver_amg(solver);
#endif

#if USE_SUPERLU
        case SLSPACK_SOLVER_SUPERLU:
            return slspack_solver_superlu(solver);
#endif

#if USE_PARDISO
        case SLSPACK_SOLVER_PARDISO:
            return slspack_solver_pardiso(solver);
#endif

#if USE_HSL_MI20
        case SLSPACK_SOLVER_MI20AMG:
            return slspack_solver_mi20amg(solver);
#endif

#if USE_SXAMG
        case SLSPACK_SOLVER_SXAMG:
            return slspack_solver_sxamg(solver);
#endif

        default:
            slspack_error(1, "solver: wrong solver type!\n");
            return -1;
    }
}

void slspack_solver_set_rtol(SLSPACK_SOLVER *s, double tol)
{
    assert(s != NULL);

    if (tol < 0) {
        slspack_warning("solver: tol is less than zero!\n");
    }
    else {
        s->tol_rel = tol;
    }
}

void slspack_solver_set_atol(SLSPACK_SOLVER *s, double tol)
{
    assert(s != NULL);

    if (tol < 0) {
        slspack_warning("solver: tol is less than zero!\n");
    }
    else {
        s->tol_abs = tol;
    }
}

void slspack_solver_set_rbtol(SLSPACK_SOLVER *s, double tol)
{
    assert(s != NULL);

    if (tol < 0) {
        slspack_warning("solver: tol is less than zero!\n");
    }
    else {
        s->tol_rbn = tol;
    }
}

void slspack_solver_set_maxit(SLSPACK_SOLVER *s, int maxit)
{
    assert(s != NULL);

    if (maxit <= 0) {
        slspack_warning("solver: maxit is less or equal to zero!\n");
    }
    else {
        s->maxit = maxit;
    }
}

void slspack_solver_set_restart(SLSPACK_SOLVER *s, int m)
{
    assert(s != NULL);

    if (m <= 0) {
        slspack_warning("solver: restart is too small!\n");
    }
    else {
        s->restart = m;
    }
}

void slspack_solver_set_verbosity(SLSPACK_SOLVER *s, int v)
{
    assert(s != NULL);

    s->verb = v;
}

double slspack_solver_get_residual(SLSPACK_SOLVER s)
{
    return s.residual;
}

int slspack_solver_get_nits(SLSPACK_SOLVER s)
{
    return s.nits;
}
