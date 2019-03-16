
#ifndef SLSPACK_DATATYPES_H
#define SLSPACK_DATATYPES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "config.h"

typedef struct SLSPACK_MAT_
{
    double *Ax;
    int *Ap;
    int *Aj;

    int num_rows;
    int num_cols;
    int num_nnzs;

} SLSPACK_MAT;

typedef struct SLSPACK_VEC_
{
    double *d;

    int n;

} SLSPACK_VEC;


/* solver type */
typedef enum SLSPACK_SOLVER_TYPE_
{
#if USE_LASPACK
    SLSPACK_SOLVER_LASPACK,   /* laspack */
#endif

#if USE_SSPARSE
    SLSPACK_SOLVER_UMFPACK,   /* ssparse, umf */
    SLSPACK_SOLVER_KLU,       /* ssparse, klu */
#endif

#if USE_MUMPS
    SLSPACK_SOLVER_MUMPS,     /* mumps */
#endif

#if USE_PETSC
    SLSPACK_SOLVER_PETSC,     /* petsc */
#endif

#if USE_LIS
    SLSPACK_SOLVER_LIS,       /* lis */
#endif

#if USE_FASP
    SLSPACK_SOLVER_FASP,      /* FASP */
    SLSPACK_SOLVER_AMG,       /* amg from FASP */
    SLSPACK_SOLVER_FMG,       /* fmg from FASP */
#endif

#if USE_SUPERLU
    SLSPACK_SOLVER_SUPERLU,   /* superlu */
#endif

#if USE_PARDISO
    SLSPACK_SOLVER_PARDISO,   /* pardiso */
#endif

#if USE_HSL_MI20
    SLSPACK_SOLVER_MI20AMG,   /* MI20 AMG */
#endif

#if USE_SXAMG
    SLSPACK_SOLVER_SXAMG,     /* SXAMG */
#endif

} SLSPACK_SOLVER_TYPE;

struct SLSPACK_SOLVER_;

typedef struct SLSPACK_SOLVER_
{
    /* universal settings */
    double tol_rel;                 /* relative tolerance */
    double tol_abs;                 /* absolute tolerance */
    double tol_rbn;                 /* relative b norm */
    int maxit;                      /* maximal number of iteration */
    int restart;                    /* used by gmres */

    /* fasp */
#if USE_FASP
    struct FASP_AMG_DATA_ *fasp_amg;
    struct FASP_DATA_     *fasp;
    int fasp_solver;
    int fasp_pc;
#endif

    /* lis */
#if USE_LIS
    void *lis;
    int lis_solver;
    int lis_pc;
#endif

    /* petsc */
#if USE_PETSC
    struct PTS_DATA_ *petsc;
    int petsc_solver;
    int petsc_pc;
#endif

    /* laspack */
#if USE_LASPACK
    struct LAS_DATA_ *laspack;
    int laspack_solver;
    int laspack_pc;
#endif

    /* mumps */
#if USE_MUMPS
    struct MPS_DATA_ *mumps;
#endif

    /* superlu */
#if USE_SUPERLU
    struct SUPERLU_DATA_ *superlu;
#endif

    /* suitesparse */
#if USE_SSPARSE
    struct SSP_DATA_ *ssparse;
#endif

#if USE_SXAMG
    struct SXAMG_DATA_ *sxamg;
#endif

    /* matrix */
    SLSPACK_MAT A;

    /* solver info */
    SLSPACK_SOLVER_TYPE type;
    SLSPACK_VEC rhs;
    SLSPACK_VEC x;

    /* return */
    double residual;                /* residual, return value */
    int nits;                       /* number of iteraion, return value */

    int verb;
    int assembled;

} SLSPACK_SOLVER;

#endif
