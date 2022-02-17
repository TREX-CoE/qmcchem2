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

AX_PROG_IRPF90
AM_PROG_AR


# Fortran compiler checks
AC_PROG_FC([ifort gfortran])
AC_FC_FREEFORM()
AC_FC_LINE_LENGTH(unlimited)
AC_FC_SRCEXT(f90)
export F77=$FC
export FFLAGS=$FCFLAGS


AC_CHECK_PROGS([QMCCHEM_MPIRUN], [srun mpiexec mpirun])
AC_SUBST([QMCCHEM_MPIRUN])

# Check for IB network interface
AX_CHECK_NIC()
QMCCHEM_NIC=$NIC
AC_SUBST([QMCCHEM_NIC])

# Checks for libraries.
AX_ZMQ([4.0], [], [ AC_MSG_ERROR([Please install ZeroMQ >= 4.0]) ])
AC_CHECK_LIB([z], [compress])
AC_CHECK_LIB([zmq], [zmq_send])
AC_CHECK_LIB([qmckl], [qmckl_context_create])
AC_CHECK_LIB([f77zmq], [f77_zmq_send__], [],
  AC_CHECK_LIB([f77zmq], [f77_zmq_send_], [],
    AC_CHECK_LIB([f77zmq], [f77_zmq_send], [], [], [-lzmq]),
  [-lzmq])
[-lzmq])
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

# Required by QMCkl
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT64_T
AC_TYPE_UINT64_T

AC_LANG([Fortran])

# Check CPU featurese
AC_ARG_WITH(arch, [AS_HELP_STRING([--with-arch],[enable avx, avx2 or avx512])], with_arch=$withval, with_arch=avx2)

AS_CASE([$FC],
  [*ifort*], [
    AS_CASE([$with_arch],
      [avx], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -xAVX"],
      [avx2], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -xCORE-AVX2"],
      [avx512], [ALIGN=64 ; FCFLAGS="${FCFLAGS} -xCORE-AVX512"])
], [*gfortran*], [
AS_CASE([$with_arch],
  [avx], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -march=sandybridge -ftree-vectorize"],
  [avx2], [ALIGN=32 ; FCFLAGS="${FCFLAGS} -march=haswell -ftree-vectorize"],
  [avx512], [ALIGN=64 ; FCFLAGS="${FCFLAGS} -march=skylake-avx512 -ftree-vectorize"])
])

AS_IF([test "x.$ALIGN" != "x."],
  IRPF90_FLAGS="--align=$ALIGN"
])
AC_SUBST([IRPF90_FLAGS])

AC_CONFIG_FILES([Makefile qmcchemrc])
AC_OUTPUT


echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION} ${QMCKL_DEVEL}

Prefix: '${prefix}'.

CC...................: ${CC}
CPPFLAGS.............: ${CPPFLAGS}
CFLAGS...............: ${CFLAGS}
FC...................: ${FC}
FCLAGS...............: ${FCFLAGS}
LDFLAGS:.............: ${LDFLAGS}
LIBS.................: ${LIBS}
CPU features.........: ${X86_FEATURE_CFLAGS}
IRPF90_FLAGS.........: ${IRPF90_FLAGS}

Package features:
  ${ARGS}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"
