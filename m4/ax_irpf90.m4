# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_irpf90.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_IRPF90
#
# DESCRIPTION
#
#   Check for the program 'irpf90', let script continue if exists, pops up
#   error message if not.
#
#   Besides checking existence, this macro also set these environment
#   variables upon completion:
#
#     IRPF90 = which irpf90
#
# DEPENDENCIES
#
#   AX_CHECK_GNU_MAKE
#
# LICENSE
#
#   Copyright (c) 2021 Anthony Scemama
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.


AU_ALIAS([AC_PROG_IRPF90], [AX_PROG_IRPF90])

AC_DEFUN([AX_PROG_IRPF90], [

# Requirements
AC_REQUIRE([AX_CHECK_GNU_MAKE])

AS_IF([test "x$ifGNUmake" = "x#"], [ AC_MSG_ERROR([GNU Make (gmake) is required with IRPF90]) ])

# IRPF90
AC_PATH_PROG([IRPF90], [irpf90], [nocommand])
AS_IF([test "x$IRPF90" = xnocommand], [
        AC_MSG_ERROR([irpf90 not found in $PATH]) ])

AC_FC_FREEFORM
AC_FC_LINE_LENGTH
AC_FC_MODULE_EXTENSION

])

