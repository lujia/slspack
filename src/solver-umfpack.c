
#include "solver-umfpack.h"

#if USE_SSPARSE

#include <umfpack.h>

typedef struct SSP_DATA_
{
    void   *sym;
    void   *num;
    int    *Ap;
    int    *Ai;
    double *Ad;
    int    factored;

} SSP_DATA;

static int umfpack_destroy(SSP_DATA *sspd)
{
    if (sspd == NULL) return 0;

    if (sspd->sym != NULL) umfpack_di_free_symbolic(&sspd->sym);
    if (sspd->num != NULL) umfpack_di_free_numeric(&sspd->num);

    slspack_free(sspd->Ap);
    slspack_free(sspd->Ai);
    slspack_free(sspd->Ad);
    slspack_free(sspd);

    return 0;
}

static int umfpack_factor(SSP_DATA *sspd, SLSPACK_SOLVER *solver)
{
    int i, j;
    int n = solver->A.num_rows, nnz = solver->A.Ap[solver->A.num_rows];
    int status;
    double time0;
    SLSPACK_MAT A = solver->A;
    int end;

    if (sspd->factored) return 0;
    sspd->factored = 1;

    /* assemble matrix */
    sspd->Ai = slspack_malloc(sizeof(int) * nnz);
    sspd->Ad = slspack_malloc(sizeof(double) * nnz);
    sspd->Ap = slspack_calloc(sizeof(int) * (n + 1));

    for (i = 0; i < n; i++) {
        end = A.Ap[i + 1];

        for (j = A.Ap[i]; j < end; j++) {
            sspd->Ap[A.Aj[j] + 1]++;
        }
    }

    for (i = 0; i < n; i++) sspd->Ap[i + 1] += sspd->Ap[i];

    for (i = 0; i < n; i++) {
        end = A.Ap[i + 1];

        for (j = A.Ap[i]; j < end; j++) {
            sspd->Ai[sspd->Ap[A.Aj[j]]] = i;
            sspd->Ad[sspd->Ap[A.Aj[j]]] = A.Ax[j];
            sspd->Ap[A.Aj[j]]++;
        }
    }

    assert(sspd->Ap[n] == nnz);
    for (i = n; i > 0; i--) sspd->Ap[i] = sspd->Ap[i - 1];
    sspd->Ap[0] = 0;

    /* symbolic fac */
    if (solver->verb > 0) slspack_printf("umfpack: analyse and factor matrix.\n");
    time0 = slspack_get_time();
    status = umfpack_di_symbolic(n, n, sspd->Ap, sspd->Ai, sspd->Ad,
            &sspd->sym, NULL, NULL) ;

    if (status != UMFPACK_OK) {
        slspack_error(1, "umfpack: error calling umfpack_di_symbolic(), status = %d.\n", status);
    }

    if (solver->verb > 0) {
        slspack_printf("umfpack: symbolic factorize: %lg s\n", slspack_get_time() - time0);
    }

    /* numeric fac */
    time0 = slspack_get_time();
    status = umfpack_di_numeric(sspd->Ap, sspd->Ai, sspd->Ad, sspd->sym, &sspd->num, NULL, NULL);

    if (status != UMFPACK_OK) {
        slspack_error(1, "umfpack: error calling umfpack_di_numeric(), status = %d.\n", status);
    }

    umfpack_di_free_symbolic(&sspd->sym) ;
    sspd->sym = NULL;

    if (solver->verb > 0) {
        slspack_printf("umfpack: numeric factorize: %lg s\n", slspack_get_time() - time0);
    }

    return 0;
}

void slspack_solver_umfpack_create(SLSPACK_SOLVER *s)
{
    assert(s != NULL);

    s->ssparse = slspack_calloc(sizeof(SSP_DATA));
}

void slspack_solver_umfpack_destroy(SLSPACK_SOLVER *s)
{
    SSP_DATA *sspd;

    if (s == NULL) return;

    sspd = s->ssparse;
    if (sspd == NULL) return;

    /* free mem */
    umfpack_destroy(sspd);
    s->ssparse = NULL;
}

int slspack_solver_umfpack(SLSPACK_SOLVER *solver)
{
    int status;
    double time0;
    SSP_DATA *sspd;
    double *x = solver->x.d;

    /* init */
    sspd = solver->ssparse;
    sspd->factored = 0;

    /* fac */
    umfpack_factor(sspd, solver);

    /* solve */
    time0 = slspack_get_time();
    status = umfpack_di_solve(UMFPACK_A, sspd->Ap, sspd->Ai, sspd->Ad, x,
            solver->rhs.d, sspd->num, NULL, NULL);

    if (status != UMFPACK_OK) {
        slspack_error(1, "umfpack: error calling umfpack_di_solve(), status = %d.\n", status);
    }

    if (solver->verb > 0) {
        slspack_printf("umfpack: solve: %lg s\n", slspack_get_time() - time0);
    }

    solver->residual = 0.0;
    return solver->nits = 1;
}

#endif
