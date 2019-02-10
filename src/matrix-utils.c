#include "matrix-utils.h"

void slspack_mat_destroy(slspack_mat_csr *A)
{
    if (A == NULL) return;

    slspack_free(A->Ap);
    slspack_free(A->Aj);
    slspack_free(A->Ax);

    bzero(A, sizeof(*A));
}

slspack_mat_csr slspack_mat_create(int nrows, int ncols, int *Ap, int *Aj, double *Ax)
{
    slspack_mat_csr A;

    assert(nrows > 0);
    assert(ncols > 0);
    assert(Ap != NULL);

    A.num_rows = nrows;
    A.num_cols = ncols;
    A.num_nnzs = Ap[nrows];

    A.Ap = slspack_copy_on(Ap, sizeof(int) * (nrows + 1));
    A.Aj = slspack_copy_on(Aj, sizeof(int) * Ap[nrows]);
    A.Ax = slspack_copy_on(Ax, sizeof(double) * Ap[nrows]);

    return A;
}

slspack_mat_csr slspack_mat_transpose(const slspack_mat_csr A)
{
    slspack_mat_csr T;

    assert(A.num_rows > 0 && A.num_cols > 0);

    T.num_rows = A.num_cols;
    T.num_cols = A.num_rows;
    T.num_nnzs = A.num_nnzs;

    if (A.num_nnzs <= 0) {
        return T;
    }

    assert(A.Ap != NULL && A.Aj != NULL && A.Ax != NULL);

    T.Ap = slspack_malloc(sizeof(int) * (T.num_rows + 1));
    T.Aj = slspack_malloc(sizeof(int) * A.num_nnzs);
    T.Ax = slspack_malloc(sizeof(double) * A.num_nnzs);

    {
        int *temp = slspack_malloc(sizeof(int) * A.num_cols);
        int num_nnzs = A.Ap[A.num_rows];
        int i, jj;

        for (i = 0; i < A.num_cols; i++) {
            temp[i] = 0;
        }

        for (i = 0; i < num_nnzs; i++)
            temp[A.Aj[i]]++;

        T.Ap[0] = 0;
        for (i = 0; i < A.num_cols; i++) {
            T.Ap[i + 1] = T.Ap[i] + temp[i];
            temp[i] = 0;
        }

        for (i = 0; i < A.num_rows; i++) {
            int row_start = A.Ap[i];
            int row_end = A.Ap[i + 1];

            for (jj = row_start; jj < row_end; jj++) {
                int col;
                int offset;

                col = A.Aj[jj];
                offset = temp[col] + T.Ap[col];

                T.Aj[offset] = i;
                T.Ax[offset] = A.Ax[jj];

                temp[col] += 1;
            }
        }

        slspack_free(temp);
    }

    return T;
}
