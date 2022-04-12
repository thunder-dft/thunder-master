#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([fireball], [1.0], [gufranco@mail.wvu.edu])
AC_CONFIG_SRCDIR([src/include/constants.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_FC
AC_PROG_INSTALL

# Checks for libraries.

# Checks for header files.

# Checks for command line options

dnl Example of default-enabled feature
dnl AC_ARG_ENABLE([foo],
dnl    AS_HELP_STRING([--disable-foo], [Disable feature foo]))
dnl
dnl AS_IF([test "x$enable_foo" != "xno"], [
dnl   dnl Do the stuff needed for enabling the feature
dnl ])

dnl THEORY
AC_ARG_ENABLE([THEORY],
    AS_HELP_STRING([--enable-THEORY=[ARG]], [Select THEORY (HARRIS (default), DOGS or KS)]))
THEORY=HARRIS
AS_IF([test "x$enable_THEORY" = "xDOGS"], [
  AC_DEFINE([DOGS], 1, [THEORY=DOGS enabled])
  THEORY=DOGS
])
AS_IF([test "x$enable_THEORY" = "xKS"], [
  AC_DEFINE([KS], 1, [THEORY=KS enabled])
  THEORY=KS
])
AC_SUBST([THEORY],[${THEORY}])

dnl METHOD
AC_ARG_ENABLE([METHOD],
    AS_HELP_STRING([--enable-METHOD=[ARG]], [Select METHOD (GAMMA (default), KPOINTS, MDET or KS)]))
METHOD=GAMMA
AS_IF([test "x$enable_METHOD" = "xKPOINTS"], [
  AC_DEFINE([KPOINTS], 1, [METHOD=KPOINTS enabled])
  METHOD=KPOINTS
])
AS_IF([test "x$enable_METHOD" = "xMDET"], [
  AC_DEFINE([MDET], 1, [METHOD=MDET enabled])
  AC_DEFINE([GAMMA], 1, [METHOD=MDET enabled and MDET is also GAMMA])
  METHOD=MDET
])
AS_IF([test "x$enable_METHOD" = "xKS"], [
  AC_DEFINE([KS], 1, [METHOD=KS enabled])
  METHOD=KS
])
AS_IF([test "x${METHOD}" = "xGAMMA"], [
  AC_DEFINE([GAMMA], 1, [METHOD=GAMMA enabled])
])

AC_SUBST([METHOD],[${METHOD}])

dnl XC
AC_ARG_ENABLE([XC],
    AS_HELP_STRING([--enable-XC=[ARG]], [Select XC (eXchange and Correlation) LDA (default), BLYP or LIBXC ]))
XC=LDA
AS_IF([test "x$enable_XC" = "xBLYP"], [
  AC_DEFINE([BLYP], 1, [XC=BLYP enabled])
  XC=BLYP
])
AS_IF([test "x$enable_XC" = "xLIBXC"], [
  AC_DEFINE([LIBXC], 1, [XC=LIBXC enabled])
  XC=LIBXC
])
AC_SUBST([XC],[${XC}])

dnl XCTHEORY
AC_ARG_ENABLE([XCTHEORY],
    AS_HELP_STRING([--enable-XCTHEORY=[ARG]], [Select XC Theory Level MCWEDA (default) or HORSFIELD ]))
XCTHEORY=MCWEDA
AS_IF([test "x$enable_XCTHEORY" = "xHORSFIELD"], [
  AC_DEFINE([HORSFIELD], 1, [XCTHEORY=HORSFIELD enabled])
  XCTHEORY=HORSFIELD
])
AC_SUBST([XCTHEORY],[${XCTHEORY}])

dnl MACHINE
AC_ARG_ENABLE([MACHINE],
    AS_HELP_STRING([--enable-MACHINE=[ARG]], [Select MACHINE (GENERIC (default), FIREBALL, SPRUCE, etc) ]))
MACHINE=GENERIC
AS_IF([test "x$enable_MACHINE" = "xFIREBALL"], [
  MACHINE=FIREBALL
])
AS_IF([test "x$enable_MACHINE" = "xSPRUCE"], [
  MACHINE=SPRUCE
])
AC_SUBST([MACHINE],[${MACHINE}])

dnl LINALG
AC_ARG_ENABLE([LINALG],
    AS_HELP_STRING([--enable-LINALG=[ARG]], [Select Linear Algebra (MKL (default), LAPACK)]))
LINALG=MKL
AS_IF([test "x$enable_LINALG" = "xLAPACK"], [
  LINALG=LAPACK
])
AC_SUBST([LINALG],[${LINALG}])

dnl OPENMP
AC_ARG_ENABLE([OPENMP],
    AS_HELP_STRING([--enable-OPENMP], [Enable OpenMP (Shared memory Multiprocessing)]))
OPENMP=no
AS_IF([test "x$enable_OPENMP" = "xyes"], [
  AC_DEFINE([OPENMP], 1, [OPENMP enabled])
  OPENMP=yes
])
AC_SUBST([OPENMP],[$OPENMP])

dnl MPI
AC_ARG_ENABLE([MPI],
    AS_HELP_STRING([--enable-MPI], [Enable MPI (Message Passing Interface)]))
MPI=no
AS_IF([test "x$enable_MPI" = "xyes"], [
  AC_DEFINE([MPI], 1, [MPI enabled])
  MPI=yes
])
AC_SUBST([MPI],[$MPI])

dnl DEBUG
AC_ARG_ENABLE([DEBUG],
    AS_HELP_STRING([--enable-DEBUG], [Enable DEBUGGING]))
DEBUG=no
AS_IF([test "x$enable_DEBUG" = "xyes"], [
  AC_DEFINE([DEBUG], 1, [DEBUG enabled])
  AC_SUBST([MODE],[DEBUG])
  DEBUG=yes
])
AC_SUBST([DEBUG],[$DEBUG])

dnl OPT
AC_ARG_ENABLE([OPT],
    AS_HELP_STRING([--enable-OPT], [Enable OPTIMIZATION]))
OPT=no
AS_IF([test "x$enable_OPT" = "xyes"], [
  AC_DEFINE([OPT], 1, [OPT enabled])
  AC_SUBST([MODE],[OPT])
  DEBUG=yes
])
AC_SUBST([OPT],[$OPT])

AC_SUBST([COMPILER],[$FC])


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/a.GLOBAL/Makefile
                 src/b.FUNCTIONS/Makefile
                 src/c.SYSTEM/Makefile
                 src/d.FUNCTIONS_EXTRA/Makefile
                 src/e.FDATA/Makefile
                 src/f.MPI/Makefile
                 src/g.XC_FUNCTIONALS/Makefile
                 src/h.SOLVESH/Makefile
                 src/j.ASSEMBLERS/Makefile
                 src/k.DASSEMBLERS/Makefile
                 src/l.SCF/Makefile
                 src/m.MD/Makefile
                 src/n.NAMD/Makefile
                 src/o.OUTPUT/Makefile
                 src/p.THEORY/Makefile
                 src/u.UTIL/Makefile])
AC_OUTPUT

echo \
"-------------------------------------------------

 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'.
 Compiler: '${FC} ${FCFLAGS} ${CPPFLAGS}'

 Package features:
   THEORY   : ${THEORY}
   METHOD   : ${METHOD}
   XC       : ${XC}
   XCTHEORY : ${XCTHEORY}

   OPENMP   : ${OPENMP}
   MPI      : ${MPI}
   DEBUG    : ${DEBUG}
   OPT      : ${OPT}

   MACHINE  : ${MACHINE}
   COMPILER : ${COMPILER}
   LINALG   : ${LINALG}
   MODE     : ${MODE}

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all                - build the convenience libraries
     install            - install everything
     clean              - remove all compiled files (.o, .a and .mod)
     cleanall           - clean + remove all binaries (.x)

--------------------------------------------------"