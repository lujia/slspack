
#include "vector.h"

SLSPACK_VEC slspack_vec_create(int n)
{
    SLSPACK_VEC v;

    assert(n >= 0);

    v.n = n;
    v.d = slspack_malloc(sizeof(*v.d) * n);

    /* init */
#if 1
    {
        int i;

        for (i = 0; i < n; i++) v.d[i] = 0.;
    }
#else
    bzero(v.d, n * sizeof(*v.d));
#endif

    return v;
}

void slspack_vec_destroy(SLSPACK_VEC *v)
{
    if (v == NULL) return;

    if (v->n > 0) slspack_free(v->d);

    v->n = 0;
    v->d = NULL;
}

/* set value */
void slspack_vec_set_value(SLSPACK_VEC x, double val)
{
    int i, n = x.n;

    for (i = 0; i < n; i++) {
        x.d[i] = val;
    }
}

/* set value */
void slspack_vec_set_value_by_array(SLSPACK_VEC x, double *val)
{
    int n = x.n;

    if (n <= 0) return;

    assert(val != NULL);
    memcpy(x.d, val, sizeof(*val) * n);
}

void slspack_vec_add_value_by_index(SLSPACK_VEC x, int i, double val)
{
    assert(i >= 0 && i < x.n);

    x.d[i] += val;
}

void slspack_vec_set_value_by_index(SLSPACK_VEC x, int i, double val)
{
    assert(i >= 0 && i < x.n);

    x.d[i] = val;
}

void slspack_vec_get_value(double *val, SLSPACK_VEC x)
{
    int n = x.n;

    assert(n >= 0);
    if (n > 0) assert(val != NULL);

    memcpy(val, x.d, sizeof(*val) * n);
}

double slspack_vec_get_value_by_index(SLSPACK_VEC x, int i)
{
    assert(i >= 0 && i < x.n);

    return x.d[i];
}

/* copy */
void slspack_vec_copy(SLSPACK_VEC x, const SLSPACK_VEC y)
{
    assert(x.n == y.n);
    assert(x.n >= 0);

    memcpy(x.d, y.d, x.n * sizeof(*x.d));
}

double slspack_vec_norm(SLSPACK_VEC x)
{
    int i, n = x.n;
    double s = 0.;

    assert(n >= 0);

    for (i = 0; i < n; i++) {
        s += x.d[i] * x.d[i];
    }

    return sqrt(s);
}
