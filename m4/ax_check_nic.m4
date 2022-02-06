# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_irpf90.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_NIC
#
# DESCRIPTION
#
#   Checks for the fastest network interface to use. Default is InfiniBand
#   (ib*) and then any other than the loopback (lo) interface. The result
#   is stored in the NIC environment variable.
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


AU_ALIAS([AC_CHECK_NIC], [AX_CHECK_NIC])


AC_DEFUN([AX_CHECK_NIC], [

  if test "$NIC.x" == ".x"; then

    # Check if ip command exists
    AC_CHECK_PROG([IP_PROG], [ip], [ip], [nocommand])

    if test "$IP_PROG.x" == "nocommand.x" ; then
      AC_MSG_ERROR([Program ip not found: Unable to detect NIC])
    fi

    INTERFACES=[`ip link show | grep -e '^[0-9]' | cut -d ':' -f 2`]
    echo $INTERFACES

    NIC=""
    for interface in $INTERFACES
    do
      case $interface in
        ib*)
          NIC=$interface
          break
          ;;
      esac
    done

    if test "$NIC.x" == ".x" ; then
      for interface in $INTERFACES
      do
        case $interface in
          lo*)
            ;;
          *)
            NIC=$interface
            break
            ;;
        esac
      done
    fi

  fi

])

