
#include "slspack.h"

static SLSPACK_MAT laplacian_5pt(const int N)
{
    SLSPACK_MAT A;
    int nz = 0;
    int i, j;

    assert (N > 0);

    A.num_rows = N * N;
    A.num_cols = N * N;
    A.num_nnzs = 5 * N * N - 4 * N;
    A.Ap = slspack_malloc(sizeof(int) * (A.num_rows + 1));
    A.Aj = slspack_malloc(sizeof(int) * A.num_nnzs);
    A.Ax = slspack_malloc(sizeof(double) * A.num_nnzs);

    A.Ap[0] = 0;
    for(i = 0; i < A.num_rows; i++) A.Ap[i + 1] = 0;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            int indx = N * i + j;

            if (i > 0){
                A.Aj[nz] = indx - N;
                A.Ax[nz] = -1;
                nz++;
            }

            if (j > 0){
                A.Aj[nz] = indx - 1;
                A.Ax[nz] = -1;
                nz++;
            }

            A.Aj[nz] = indx;
            A.Ax[nz] = 4;
            nz++;

            if (j < N - 1){
                A.Aj[nz] = indx + 1;
                A.Ax[nz] = -1;
                nz++;
            }

            if (i < N - 1){
                A.Aj[nz] = indx + N;
                A.Ax[nz] = -1;
                nz++;
            }

            A.Ap[indx + 1] = nz;
        }
    }

    return A;
}

int main(void)
{
    int i, n;
    int dim = 100;
    int m = 30;
    int itr_max = 3000;

    SLSPACK_MAT A;
    SLSPACK_SOLVER solver;
    slspack_vec xg;
    slspack_vec bc;
    double tc, ta;

    slspack_printf("CSR: laplacian, grid size %d\n", dim);
    A = laplacian_5pt(dim);

    n = A.num_rows;
    dim = A.num_nnzs;
    slspack_printf("CSR: rows: %d", n);
    slspack_printf(" nonzeros: %d mem (csr): %.3f Mb\n", dim, \
            ((dim + n) * sizeof(int) + dim * sizeof(double))/ 1048576.);

    /* generate vector x, b */
    xg = slspack_vec_create(n);
    bc = slspack_vec_create(n);

    for (i = 0; i < n; i++) {
        slspack_vec_set_value_by_index(xg, i, 0.);
        slspack_vec_set_value_by_index(bc, i, 1.);
    }

    /* step 1: set up solver */
    slspack_solver_create(&solver, SLSPACK_SOLVER_AMG);

    /* step 2: change default settings (optional) */
    slspack_solver_set_restart(&solver, m);
    slspack_solver_set_maxit(&solver, itr_max);
    slspack_solver_set_verbosity(&solver, 8);

    /* step 3: assemble solver and preconditioner */
    slspack_solver_assemble(&solver, A, xg, bc);

    /* step 4: solve */
    ta = slspack_get_time();
    slspack_solver_solve(&solver);
    tc = slspack_get_time() - ta;

    slspack_printf("total solver time: %g \n", tc);
    slspack_printf("solution L2 norm: %g residual: %g\n", slspack_vec_norm(xg), solver.residual);

    /* step 5: destroy solver */
    slspack_solver_destroy(&solver);

    slspack_mat_destroy(&A);
    slspack_vec_destroy(&xg);
    slspack_vec_destroy(&bc);

    return 0;
}
