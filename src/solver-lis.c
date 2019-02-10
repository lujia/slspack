
#include "solver-lis.h"

#if USE_LIS

#include "lis.h"

static char * lis_solver[] = {
    "-i bicgstab",
    "-i bicgstabl",
    "-i cg",
    "-i cgs",
    "-i bicg",
    "-i bicgsafe",
    "-i bicr",
    "-i cr",
    "-i bicrstab",
    "-i bicrsafe",
    "-i idrs",
    "-i crs",
    "-i gpbicr",
    "-i gpbicg",
    "-i tfqmr",
    "-i orthomin",
    "-i gmres",
    "-i fgmres",
    "-i minres",
};

static char * lis_pc[] = {
    "-p none",
    "-p ilut",
    "-p ilu -ilu_fill 1",
    "-p is",
    "-p sainv",
    "-p saamg -saamg_unsym -saamg_theta 0.5",
    "-p hybrid",
    "-p iluc",
    "-p ssor",
    "-p jacobi",
};

#define Nlis_solver (sizeof(lis_solver)/sizeof(lis_solver[0]))
#define Nlis_pc     (sizeof(lis_pc)/sizeof(lis_pc[0]))

void slspack_solver_lis_create(SLSPACK_SOLVER *solver)
{
    assert(solver != NULL);
    solver->lis = slspack_malloc(sizeof(LIS_SOLVER));

    lis_solver_create(solver->lis);
}

void slspack_solver_lis_destroy(SLSPACK_SOLVER *solver)
{
    LIS_SOLVER s;

    if (solver == NULL || solver->lis == NULL) return;

    s = *(LIS_SOLVER *)solver->lis;
    lis_solver_destroy(s);

    slspack_free(solver->lis);
    solver->lis = NULL;
}

int slspack_solver_lis(SLSPACK_SOLVER *solver)
{
    double time0;
    double *x = solver->x.d;
    slspack_mat_csr A = solver->A;

    LIS_MATRIX lA;
    LIS_VECTOR lx, lb;
    LIS_SCALAR t;
    LIS_INT iter;
    int lis_solver_id = solver->lis_solver;
    int lis_pc_id = solver->lis_pc;

    /* init */
    lis_matrix_create(LIS_COMM_WORLD,&lA);
    lis_matrix_set_size(lA, 0, A.num_rows);
    lis_vector_duplicate(lA, &lx);
    lis_vector_duplicate(lA, &lb);

    /* assemble */
    for (int i = 0; i < A.num_rows; i++) {
        int end = A.Ap[i + 1];

        for (int j = A.Ap[i]; j < end; j++) {
            lis_matrix_set_value(LIS_INS_VALUE, i, A.Aj[j], A.Ax[j], lA);
        }

        lis_vector_set_value(LIS_INS_VALUE, i, x[i], lx);
        lis_vector_set_value(LIS_INS_VALUE, i, solver->rhs.d[i], lb);
    }

    lis_matrix_assemble(lA);

    /* solve */
    time0 = slspack_get_time();

    /* check */
    if ((int)Nlis_solver <= lis_solver_id || lis_solver_id < 0) {
        slspack_warning("lis: wrong solver id: %d, set to 0\n", lis_solver_id);
        lis_solver_id = 0;
    }

    if ((int)Nlis_pc <= lis_pc_id || lis_pc_id < 0) {
        slspack_warning("lis: wrong pc id: %d, set to 2\n", lis_pc_id);
        lis_pc_id = 2;
    }

    {
        char sn[256];

        sprintf(sn, "-tol %e", solver->tol_rel);
        lis_solver_set_option(sn, *(LIS_SOLVER *)solver->lis);

        sprintf(sn, "-maxiter %d", solver->maxit);
        lis_solver_set_option(sn, *(LIS_SOLVER *)solver->lis);

        if (solver->verb > 1) {
            lis_solver_set_option("-print 2", *(LIS_SOLVER *)solver->lis);
        }
    }

    lis_solver_set_option(lis_solver[lis_solver_id], *(LIS_SOLVER *)solver->lis);
    lis_solver_set_option(lis_pc[lis_pc_id], *(LIS_SOLVER *)solver->lis);

    lis_solve(lA, lb, lx, *(LIS_SOLVER *)solver->lis);

    if (solver->verb > 2) {
        slspack_printf("lis: solve: %lg s\n", slspack_get_time() - time0);
    }

    /* return */
    for (int i = 0; i < A.num_rows; i++) {
        lis_vector_get_value(lx, i, &t);
        x[i] = t;
    }

    lis_solver_get_iter(*(LIS_SOLVER *)solver->lis, &iter);
    lis_solver_get_residualnorm(*(LIS_SOLVER *)solver->lis, &t);

    /* free mem */
    lis_matrix_destroy(lA);
    lis_vector_destroy(lx);
    lis_vector_destroy(lb);

    /* relative error */
    solver->residual = t;
    solver->nits = iter;

    return iter;
}

void slspack_solver_lis_set_pars(SLSPACK_SOLVER *solver, unsigned int solver_id, unsigned int pc_id)
{
    assert(solver != NULL);

    solver->lis_solver = solver_id;
    solver->lis_pc = pc_id;
}

void slspack_solver_lis_set_option(SLSPACK_SOLVER *solver, char *o)
{
    assert(solver != NULL);
    if (solver->lis == NULL) return;

    lis_solver_set_option(o, *(LIS_SOLVER *)solver->lis);
}

#endif
