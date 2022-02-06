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

AC_CHECK_PROGS([QMCCHEM_MPIRUN], [srun mpiexec mpirun])
AC_SUBST([QMCCHEM_MPIRUN])

# Check for IB network interface
AX_CHECK_NIC()
QMCCHEM_NIC=$NIC
echo "NIC=$NIC"
echo "QMCCHEM_NIC=$QMCCHEM_NIC"
AC_SUBST([QMCCHEM_NIC])

# Fortran compiler checks
AC_PROG_FC([ifort gfortran])
AC_FC_LINE_LENGTH([132])


# Checks for libraries.
AX_ZMQ([4.0], [], [ AC_MSG_ERROR([Please install ZeroMQ >= 4.0]) ])
AX_BLAS()
AX_LAPACK()

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
#AC_CHECK_HEADER_STDBOOL
#AC_TYPE_INT64_T
#AC_TYPE_UINT64_T


AC_LANG([Fortran])

AC_CONFIG_FILES([Makefile qmcchemrc])
AC_OUTPUT
