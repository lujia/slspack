
#include "solver-amg.h"

#if USE_FASP

typedef struct FASP_AMG_DATA_
{
    AMG_param pars;

} FASP_AMG_DATA;

void slspack_solver_amg_create(SLSPACK_SOLVER *s)
{
    assert(s != NULL);

    s->fasp_amg = slspack_malloc(sizeof(FASP_AMG_DATA));
    fasp_param_amg_init(&s->fasp_amg->pars);
}

void slspack_solver_amg_destroy(SLSPACK_SOLVER *s)
{
    if (s == NULL) return;

    if (s->fasp_amg != NULL) {
        slspack_free(s->fasp_amg);
        s->fasp_amg = NULL;
    }
}

/* algebraic multi-grid solver */
static void slspack_solver_amg_assemble(dCSRmat *A, dvector *b, dvector *x, SLSPACK_SOLVER *solver)
{
    int i, n;
    INT nz;
    
    /* mat */
    n = solver->A.num_rows;
    A->row = n;
    A->col = n;
    A->IA  = (INT *)fasp_mem_calloc(n+1, sizeof(INT));
    
    for ( i = 0; i <= n; ++i ) {
        A->IA[i] = solver->A.Ap[i];
    }
    
    nz = A->IA[n];
    A->nnz = nz;
    A->JA  = (INT *) fasp_mem_calloc(nz, sizeof(INT));
    A->val = (REAL *)fasp_mem_calloc(nz, sizeof(REAL));
    
    for ( i = 0; i < nz; ++i ) {
        A->JA[i] = solver->A.Aj[i];
        A->val[i] = solver->A.Ax[i];
    }
    
    /* b and x */
    b->row = n;
    b->val = (REAL *)fasp_mem_calloc(n, sizeof(REAL));

    x->row = n;
    x->val = (REAL *)fasp_mem_calloc(n, sizeof(REAL));
    
    for ( i = 0; i < n; ++i ) {
        b->val[i] = solver->rhs.d[i];
        x->val[i] = solver->x.d[i];
    }
}

int slspack_solver_amg(SLSPACK_SOLVER *solver)
{
    dCSRmat A;
    dvector b, x;
    int status = FASP_SUCCESS;
    int i;

    int print_level;
    double time = slspack_get_time();
    AMG_param solver_amg_pars_ = solver->fasp_amg->pars;

    // set solver parameters
    solver_amg_pars_.maxit = solver->maxit;
    solver_amg_pars_.tol = solver->tol_rel;

    print_level = PRINT_SOME;
    if (solver->verb <= 0) {
        print_level = PRINT_NONE;
    }
    else if (solver->verb <= 1) {
        print_level = PRINT_SOME;
    }
    else if (solver->verb <= 2) {
        print_level = PRINT_MORE;
    }
    else if (solver->verb <= 3) {
        print_level = PRINT_MOST;
    }
    else {
        print_level = PRINT_ALL;
    }

    solver_amg_pars_.print_level = print_level;

    // read A and b 
    slspack_solver_amg_assemble(&A, &b, &x, solver);

    // solve the system //
    if (solver->type == SLSPACK_SOLVER_AMG) {      // AMG as the iterative solver
        if (print_level > PRINT_NONE) fasp_param_amg_print(&solver_amg_pars_);

        fasp_solver_amg(&A, &b, &x, &solver_amg_pars_);
    }
    else if (solver->type == SLSPACK_SOLVER_FMG) { // Full AMG as the iterative solver 
        if (print_level > PRINT_NONE) fasp_param_amg_print(&solver_amg_pars_);
        fasp_solver_famg(&A, &b, &x, &solver_amg_pars_);
    }
    else {
        slspack_printf("slspack: fasp, wrong solver type!!!\n");
        status = ERROR_SOLVER_TYPE;
        return status;
    }

    if (status < 0) {
        slspack_printf("slspack: fasp, solver failed! Exit status = %d.\n\n", status);
    }

    for ( i = 0; i < b.row; ++i ) solver->x.d[i] = x.val[i];

    time = slspack_get_time() - time;
    if (solver->verb >= 2) {
        slspack_printf("slspack: fasp, total time: %g\n", time);
    }

    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&x);

    return status;
}

void slspack_solver_amg_set_pars(SLSPACK_SOLVER *s, AMG_param par)
{
    assert(s != NULL);

    if (s->fasp_amg != NULL) {
        s->fasp_amg->pars = par;
    }
}
#endif
