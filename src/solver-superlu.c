
#include "solver-superlu.h"

#if USE_SUPERLU

#include <slu_ddefs.h>

typedef struct SUPERLU_DATA_
{
    superlu_options_t options;

} SUPERLU_DATA;

void slspack_solver_superlu_create(SLSPACK_SOLVER *s)
{
    assert(s != NULL);

    s->superlu = slspack_malloc(sizeof(SUPERLU_DATA));
    set_default_options(&s->superlu->options);
}

void slspack_solver_superlu_destroy(SLSPACK_SOLVER *s)
{
    if (s == NULL || s->superlu == NULL) return;

    slspack_free(s->superlu);
    s->superlu = NULL;
}

int slspack_solver_superlu(SLSPACK_SOLVER *solver)
{
    SuperMatrix A, L, U, B;
    slspack_mat_csr Ai = solver->A;
    int *perm_r;
    int *perm_c, i;
    int info;
    int m = Ai.num_rows, n = Ai.num_cols, nnz = Ai.num_nnzs;
    SuperLUStat_t stat;
    double time0 = slspack_get_time();
    double *BB;
    int *Aj, *Ap;
    double *b, *Ax;

    /* copy */
    Aj = slspack_copy_on(Ai.Aj, sizeof(int) * nnz);
    Ap = slspack_copy_on(Ai.Ap, sizeof(int) * (n + 1));
    Ax = slspack_copy_on(Ai.Ax, sizeof(double) * nnz);
    b = slspack_copy_on(solver->rhs.d, sizeof(double) * n);

    /* matrix A */
    assert(m == n);
    dCreate_CompCol_Matrix(&A, m, n, nnz, Ax, Aj, Ap, SLU_NR, SLU_D, SLU_GE);

    /* right-hand side */
    dCreate_Dense_Matrix(&B, m, 1, b, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    solver->superlu->options.ColPerm = COLAMD;
    StatInit(&stat);

    /* SuperLU */
    dgssv(&solver->superlu->options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    BB = (double *)((DNformat *) B.Store)->nzval;
    for (i = 0; i < n; i++) solver->x.d[i] = BB[i];

    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
    slspack_free(b);

    if (solver->verb > 0) {
        slspack_printf("superlu: solve: %lg s\n", slspack_get_time() - time0);
    }

    solver->residual = 0.0;
    return 1; 

}

#endif
