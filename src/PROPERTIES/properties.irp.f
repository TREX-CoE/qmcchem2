BEGIN_SHELL [ /usr/bin/env python3 ]
import os
from properties import properties
root = os.environ['QMCCHEM_PATH']

template = """
BEGIN_PROVIDER [ logical, calc_%(p)s ]
  implicit none
  BEGIN_DOC
  ! If true, calculate %(p)s
  END_DOC
  calc_%(p)s = .False.
  logical :: has_%(p)s
  if (.not.is_worker) then
    call ezfio_has_properties_%(p)s(has_%(p)s)
    if (has_%(p)s) then
      call ezfio_get_properties_%(p)s(calc_%(p)s)
    endif
  else
    call zmq_ezfio_has('properties_%(p)s',has_%(p)s)
    if (has_%(p)s) then
      call zmq_ezfio_get_logical('properties_%(p)s',calc_%(p)s,1)
    endif
  endif
END_PROVIDER
"""

for p in properties:
  print (template%{'p':p[1]})

t="""
 BEGIN_PROVIDER [ $T, $X_min ]
&BEGIN_PROVIDER [ $T, $X_max ]
  implicit none
  BEGIN_DOC
  ! Minimum and maximum values of $X
  END_DOC
  $X_min = huge(1.)
  $X_max =-huge(1.)
END_PROVIDER

BEGIN_PROVIDER [ $T, $X_2 $D ]
  implicit none
  BEGIN_DOC
  ! Square of $X
  END_DOC
  $X_2= $X*$X
END_PROVIDER
"""

for p in properties:
  d = ""
  if p[2] != '':
    d = ", %s"%(p[2])
  print(t.replace("$T",p[0]).replace("$X",p[1]).replace("$D",d))
END_SHELL

!==========================================================================!
! DIMENSIONS
!==========================================================================!

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
make_dims()
END_SHELL

!==========================================================================!
!                                                                         !
!==========================================================================!


!==========================================================================!
! PROPERTIES                                                              !
!==========================================================================!



