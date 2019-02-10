
#include "vector.h"

slspack_vec slspack_vec_create(int n)
{
    slspack_vec v;

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

void slspack_vec_destroy(slspack_vec *v)
{
    if (v == NULL) return;

    if (v->n > 0) slspack_free(v->d);

    v->n = 0;
    v->d = NULL;
}

/* set value */
void slspack_vec_set_value(slspack_vec x, double val)
{
    int i, n = x.n;

    for (i = 0; i < n; i++) {
        x.d[i] = val;
    }
}

/* set value */
void slspack_vec_set_value_by_array(slspack_vec x, double *val)
{
    int n = x.n;

    if (n <= 0) return;

    assert(val != NULL);
    memcpy(x.d, val, sizeof(*val) * n);
}

void slspack_vec_add_value_by_index(slspack_vec x, int i, double val)
{
    assert(i >= 0 && i < x.n);

    x.d[i] += val;
}

void slspack_vec_set_value_by_index(slspack_vec x, int i, double val)
{
    assert(i >= 0 && i < x.n);

    x.d[i] = val;
}

void slspack_vec_get_value(double *val, slspack_vec x)
{
    int n = x.n;

    assert(n >= 0);
    if (n > 0) assert(val != NULL);

    memcpy(val, x.d, sizeof(*val) * n);
}

double slspack_vec_get_value_by_index(slspack_vec x, int i)
{
    assert(i >= 0 && i < x.n);

    return x.d[i];
}

/* copy */
void slspack_vec_copy(slspack_vec x, const slspack_vec y)
{
    assert(x.n == y.n);
    assert(x.n >= 0);

    memcpy(x.d, y.d, x.n * sizeof(*x.d));
}

double slspack_vec_norm(slspack_vec x)
{
    int i, n = x.n;
    double s = 0.;

    assert(n >= 0);

    for (i = 0; i < n; i++) {
        s += x.d[i] * x.d[i];
    }

    return sqrt(s);
}
