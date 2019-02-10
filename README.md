# SLSPack

SLSPack is an interface library for famous Sparse Linear Solver PACKages, 
which has implemented interfaces for
over ten solver packages. By using this library, users can access more than ten packages, including
iterative solvers and direct solvers. This package does not implement any internal solver or preconditioner.

# Credit
All solver packages are developed by other people, such as Lis, FASP, SuperLU, PETSc, MUMPS and PARDISO.
Please cite them properly if you use those packages.

# Purpose
The purpose is to track the progress of other solver packages.

# How-to
## Configuration
The simplest way to configure is to run command:
```
./configure
```
This command will try to find optional packages from certain directories. Searching details are included by configure.in and some are explained below.

## Options
Run command:
```
./configure --help
```

Output:
```
`configure' configures this package to adapt to many kinds of systems.

Usage: ./configure [OPTION]... [VAR=VALUE]...

To assign environment variables (e.g., CC, CFLAGS...), specify them as
VAR=VALUE.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

Configuration:
  -h, --help              display this help and exit
      --help=short        display options specific to this package
      --help=recursive    display the short help of all the included packages
  -V, --version           display version information and exit
  -q, --quiet, --silent   do not print `checking ...' messages
      --cache-file=FILE   cache test results in FILE [disabled]
  -C, --config-cache      alias for `--cache-file=config.cache'
  -n, --no-create         do not create output files
      --srcdir=DIR        find the sources in DIR [configure dir or `..']

Installation directories:
  --prefix=PREFIX         install architecture-independent files in PREFIX
                          [/usr/local/slspack]
  --exec-prefix=EPREFIX   install architecture-dependent files in EPREFIX
                          [PREFIX]

By default, `make install' will install all the files in
`/usr/local/slspack/bin', `/usr/local/slspack/lib' etc.  You can specify
an installation prefix other than `/usr/local/slspack' using `--prefix',
for instance `--prefix=$HOME'.

For better control, use the options below.

Fine tuning of the installation directories:
  --bindir=DIR            user executables [EPREFIX/bin]
  --sbindir=DIR           system admin executables [EPREFIX/sbin]
  --libexecdir=DIR        program executables [EPREFIX/libexec]
  --sysconfdir=DIR        read-only single-machine data [PREFIX/etc]
  --sharedstatedir=DIR    modifiable architecture-independent data [PREFIX/com]
  --localstatedir=DIR     modifiable single-machine data [PREFIX/var]
  --libdir=DIR            object code libraries [EPREFIX/lib]
  --includedir=DIR        C header files [PREFIX/include]
  --oldincludedir=DIR     C header files for non-gcc [/usr/include]
  --datarootdir=DIR       read-only arch.-independent data root [PREFIX/share]
  --datadir=DIR           read-only architecture-independent data [DATAROOTDIR]
  --infodir=DIR           info documentation [DATAROOTDIR/info]
  --localedir=DIR         locale-dependent data [DATAROOTDIR/locale]
  --mandir=DIR            man documentation [DATAROOTDIR/man]
  --docdir=DIR            documentation root [DATAROOTDIR/doc/PACKAGE]
  --htmldir=DIR           html documentation [DOCDIR]
  --dvidir=DIR            dvi documentation [DOCDIR]
  --pdfdir=DIR            pdf documentation [DOCDIR]
  --psdir=DIR             ps documentation [DOCDIR]

System types:
  --build=BUILD     configure for building on BUILD [guessed]
  --host=HOST       cross-compile to build programs to run on HOST [BUILD]

Optional Features:
  --disable-option-checking  ignore unrecognized --enable/--with options
  --disable-FEATURE       do not include FEATURE (same as --enable-FEATURE=no)
  --enable-FEATURE[=ARG]  include FEATURE [ARG=yes]
  --enable-rpath          enable use of rpath (default)
  --disable-rpath         disable use of rpath
  --with-rpath-flag=FLAG  compiler flag for rpath (e.g., "-Wl,-rpath,")
  --enable-blas         enable BLAS support (default)
  --disable-blas        disable BLAS support
  --with-blas=blas BLAS lib
  --enable-lapack         enable LAPACK support (default)
  --disable-lapack        disable LAPACK support
  --with-lapack=lapack LAPACK lib
  --enable-laspack      enable LASPACK support (default)
  --disable-laspack     disable LASPACK support
  --with-laspack-libdir=DIR path for LASPACK library
  --with-laspack-incdir=DIR path for LASPACK header file
  --enable-ssparse      enable SSPARSE support (default)
  --disable-ssparse     disable SSPARSE support
  --with-ssparse-libdir=DIR path for SSPARSE library
  --with-ssparse-incdir=DIR path for SSPARSE header file
  --enable-petsc        enable PETSC solver (default)
  --disable-petsc       disable PETSC solver
  --with-petsc-incdir=DIR PETSC header files directory
  --with-petsc-libdir=DIR PETSC libraries directory
  --enable-mumps        enable MUMPS solver (default)
  --disable-mumps       disable MUMPS solver
  --with-mumps-incdir=DIR MUMPS header files directory
  --with-mumps-libdir=DIR MUMPS libraries directory
  --enable-lis      enable LIS support (default)
  --disable-lis     disable LIS support
  --with-lis-libdir=DIR path for LIS library
  --with-lis-incdir=DIR path for LIS header file
  --enable-fasp      enable FASP support (default)
  --disable-fasp     disable FASP support
  --with-fasp-libdir=DIR path for FASP library
  --with-fasp-incdir=DIR path for FASP header file
  --enable-superlu      enable SUPERLU support (default)
  --disable-superlu     disable SUPERLU support
  --with-superlu-libdir=DIR path for SUPERLU library
  --with-superlu-incdir=DIR path for SUPERLU header file
  --enable-pardiso      enable PARDISO support (default)
  --disable-pardiso     disable PARDISO support
  --with-pardiso-libdir=DIR path for PARDISO library
  --with-pardiso-incdir=DIR path for PARDISO header file
  --enable-hslmi20      enable HSL_MI20 support (default)
  --disable-hslmi20     disable HSL_MI20 support
  --with-hslmi20-libdir=DIR path for HSL_MI20 library
  --with-hslmi20-incdir=DIR path for HSL_MI20 header file
  --enable-sxamg      enable SXAMG support (default)
  --disable-sxamg     disable SXAMG support
  --with-sxamg-libdir=DIR path for SXAMG library
  --with-sxamg-incdir=DIR path for SXAMG header file

Some influential environment variables:
  CC          C compiler command
  CFLAGS      C compiler flags
  LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
              nonstandard directory <lib dir>
  LIBS        libraries to pass to the linker, e.g. -l<library>
  CPPFLAGS    (Objective) C/C++ preprocessor flags, e.g. -I<include dir> if
              you have headers in a nonstandard directory <include dir>
  FC          Fortran compiler command
  FCFLAGS     Fortran compiler flags
  CPP         C preprocessor

Use these variables to override the choices made by `configure' or to help
it to find libraries and programs with nonstandard names/locations.

Report bugs to the package provider.
```

Most optional packages are enabled by default. However, a package can be disabled when configuring, such as "--disable-fasp" to disable FASP. When a package needs an include path and a library path, they can be set by configure. If configure cannot find correct paths, users can help configure by using options.


## Compilation
After configuration, a simple **make** command can compile the package:
```
make
```

## Installation
Run command:
```
make install
```
The package will be installed to a directory. The default is /usr/local/slspack/. A different directory can be set by --prefix=DIR.

# Solver Types

```
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

#if USE_SX_AMG
    SLSPACK_SOLVER_SXAMG,     /* SXAMG */
#endif

} SLSPACK_SOLVER_TYPE;
```
