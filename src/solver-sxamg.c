
#include "solver-sxamg.h"

#if USE_SXAMG

typedef struct SXAMG_DATA_
{
    SX_AMG_PARS pars;

} SXAMG_DATA;

void slspack_solver_sxamg_create(SLSPACK_SOLVER *s)
{
    assert(s != NULL);

    s->sxamg = slspack_malloc(sizeof(SXAMG_DATA));
    sx_amg_pars_init(&s->sxamg->pars);
}

void slspack_solver_sxamg_destroy(SLSPACK_SOLVER *s)
{
    if (s != NULL && s->sxamg != NULL) {
        slspack_free(s->sxamg);
        s->sxamg = NULL;
    }
}

int slspack_solver_sxamg(SLSPACK_SOLVER *solver)
{
    SLSPACK_MAT A = solver->A;
    SX_INT *Aj, *Ap, nnz, n, i;
    SX_FLOAT *Ax;
    SX_MAT XA;
    SX_VEC b, x;
    SX_RTN rtn;

    /* matrix */
    n = A.num_rows;
    nnz = A.num_nnzs;

    if (sizeof(int) == sizeof(SX_INT)) {
        Ap = (SX_INT *) A.Ap;
        Aj = (SX_INT *) A.Aj;
    }
    else {
        Ap = slspack_malloc(sizeof(SX_INT) * (n + 1));
        Aj = slspack_malloc(sizeof(SX_INT) * nnz);

        for (i = 0; i <= A.num_rows; i++) Ap[i] = A.Ap[i];

        for (i = 0; i < A.Ap[n]; i++) Aj[i] = A.Aj[i];
    }

    if (sizeof(double) == sizeof(SX_FLOAT)) {
        Ax = (SX_FLOAT *) A.Ax;
    }
    else {
        Ax = slspack_malloc(sizeof(SX_FLOAT) * nnz);
        for (i = 0; i < A.Ap[n]; i++) Ax[i] = A.Ax[i];
    }

    XA = sx_mat_create(n, n, Ap, Aj, Ax);

    /* vector */
    b = sx_vec_create(n);
    x = sx_vec_create(n);

    for (i = 0; i < n; i++) {
        sx_vec_set_entry(&x, i, solver->x.d[i]);
        sx_vec_set_entry(&b, i, solver->rhs.d[i]);
    }

    /* pars */
    solver->sxamg->pars.maxit = solver->maxit;
    solver->sxamg->pars.verb = solver->verb;
    solver->sxamg->pars.tol = solver->tol_rel;

    /* solve */
    rtn = sx_solver_amg(&XA, &x, &b, &solver->sxamg->pars);

    /* return solution */
    for (i = 0; i < n; i++) {
        solver->x.d[i] = sx_vec_get_entry(&x, i);
    }

    sx_mat_destroy(&XA);
    sx_vec_destroy(&b);
    sx_vec_destroy(&x);

    if (sizeof(int) != sizeof(SX_INT)) {
        slspack_free(Ap);
        slspack_free(Aj);
    }

    if (sizeof(double) != sizeof(SX_FLOAT)) {
        slspack_free(Ax);
    }

    /* abs res */
    solver->residual = rtn.ares;
    solver->nits = rtn.nits;

    return rtn.nits;
}

void slspack_solver_sxamg_set_pars(SLSPACK_SOLVER *solver, SX_AMG_PARS *pars)
{
    if (solver == NULL) return;
    if (solver->sxamg == NULL) return;

    if (pars != NULL) solver->sxamg->pars = *pars;
}

#endif
