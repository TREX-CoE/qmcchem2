#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([QMC=Chem], [2.0.0], [https://gitlab.com/scemama/qmcchem/-/issues/new])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([foreign subdir-objects silent-rules])

AC_CONFIG_SRCDIR([README.md])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([include/config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AM_PROG_AR
PKG_PROG_PKG_CONFIG

AC_CHECK_PROGS([OCAML], [ocamlc], [no])
if test "$OCAML" == "no" ; then
  AC_MSG_ERROR([OCaml compiler not found])
fi

AC_CHECK_PROGS([QP], [qp_run], [no])
if test "$QP" == "no" ; then
  AC_MSG_ERROR([Quantum Package is not installed or not sourced in the current shell])
fi


AC_CHECK_PROGS([OCAMLBUILD], [ocamlbuild], [no])
if test "$OCAMLBUILD" == "no"; then
 AC_MSG_ERROR([ocamlbuild not found])
fi


# Fortran compiler checks
AC_PROG_FC([ifort gfortran])
AC_FC_FREEFORM()
AC_FC_LINE_LENGTH(unlimited)
AC_FC_SRCEXT(f90)
export F77=$FC
export FFLAGS=$FCFLAGS

AS_IF([test x"$IRPF90" == x""], [IRPF90=$PWD/external/irpf90/bin/irpf90])
AX_PROG_IRPF90

AC_CHECK_PROGS([QMCCHEM_MPIRUN], [srun mpiexec mpirun])
AC_SUBST([QMCCHEM_MPIRUN])

# Check for IB network interface
AX_CHECK_NIC()
QMCCHEM_NIC=$NIC
AC_SUBST([QMCCHEM_NIC])

# Checks for libraries.
AC_LANG_PUSH(C)
AX_ZMQ([], [], [AC_MSG_ERROR([ZeroMQ not found]) ])
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([z], [compress])
AC_CHECK_LIB([zmq], [zmq_send])



PKG_CHECK_EXISTS([qmckl], [
  PKG_CHECK_MODULES([QMCKL], [qmckl],
    [LIBS="$LIBS $QMCKL_LIBS"
     CFLAGS="$CFLAGS $QMCKL_CFLAGS"])
])
AC_CHECK_LIB([qmckl], [qmckl_context_create], [],
 [AC_MSG_ERROR([qmckl not found])])

PKG_CHECK_EXISTS([trexio], [
  PKG_CHECK_MODULES([TREXIO], [trexio],
    [LIBS="$LIBS $TREXIO_LIBS"
     CFLAGS="$CFLAGS $TREXIO_CFLAGS"])
])
AC_CHECK_LIB([trexio], [trexio_read_mo_num_32], [],
 [AC_MSG_ERROR([TREXIO not found])])

AC_LANG_PUSH([Fortran])
AC_CHECK_LIB([f77zmq], [f77_zmq_send], [], [AC_MSG_ERROR([f77zmq not found])])
AC_LANG_POP()

AS_IF([test $FC == ifort],
  [BLAS_LIBS="-mkl=sequential"], [
  AX_BLAS()
  AX_LAPACK()
  BLAS_LIBS="$LAPACK_LIBS $BLAS_LIBS"
])
LIBS="$LIBS $BLAS_LIBS"

# Required by EZFIO
AC_CHECK_FUNCS([mkdir strerror])
AC_TYPE_SIZE_T


# Required by OCaml C bindings
AC_CHECK_FUNCS([inet_ntoa])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([socket])
AC_CHECK_HEADERS([arpa/inet.h netinet/in.h sys/ioctl.h sys/socket.h])

# Required by ZeroMQ
AC_CHECK_HEADERS([stddef.h])
AC_TYPE_INT32_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Required by QMCkl and trexio
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_CHECK_HEADERS([qmckl_f.F90])
AC_CHECK_HEADERS([trexio_f.f90])

AC_LANG([Fortran])

# Add better compiler flags with ifort
AS_CASE([$FC],
   [*ifort*], [FCFLAGS="$FCFLAGS -ip -ftz -finline"],
])

# Check CPU features
AC_ARG_WITH(arch, [AS_HELP_STRING([--with-arch],[enable avx, avx2 or avx512])], with_arch=$withval, with_arch=unknown)
ALIGN=32

AS_CASE([$FC],
  [*ifort*], [
    AS_CASE([$with_arch],
      [avx], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -march=corei7-avx"],
      [avx2], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -march=core-avx2"],
      [avx512], [ALIGN=64 ; FCFLAGS="${FCFLAGS} -march=skylake-avx512"])
], [*gfortran*], [
AS_CASE([$with_arch],
  [avx], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -march=sandybridge -ftree-vectorize"],
  [avx2], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -march=haswell -ftree-vectorize"],
  [avx512], [ALIGN=64 ; FCFLAGS="${FCFLAGS} -march=skylake-avx512 -ftree-vectorize"])
])

AS_IF([test "x.$ALIGN" != "x."], [IRPF90_FLAGS="--align=$ALIGN"],[IRPF90_FLAGS="--align=32"])
AC_SUBST([IRPF90_FLAGS])
AC_SUBST([IRPF90])

AC_CONFIG_FILES([Makefile qmcchemrc])
AC_OUTPUT


echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'.

CC...................: ${CC}
CFLAGS...............: ${CFLAGS}
FC...................: ${FC}
FCLAGS...............: ${FCFLAGS}
LIBS.................: ${LIBS}
LDFLAGS..............: ${LDFLAGS}
CPU features.........: ${with_arch}
IRPF90...............: ${IRPF90}
IRPF90_FLAGS.........: ${IRPF90_FLAGS}
OCAML................: ${OCAML}
OCAMLBUILD...........: ${OCAMLBUILD}

Package features:
  ${ARGS}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"

# With Intel, default flags should be set to
# -O2 -g -ip -ftz -finline 
