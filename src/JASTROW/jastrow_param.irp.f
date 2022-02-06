! Input data
! ----------

BEGIN_PROVIDER  [ logical, do_jast ]
 implicit none
 BEGIN_DOC  
! If true, compute the Jastrow factor
 END_DOC
 include '../types.F'
 do_jast = jast_type /= t_None
 call linfo(irp_here,'do_jast',do_jast)
END_PROVIDER

BEGIN_PROVIDER [ integer, jast_type ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Type of Jastrow factor : Simple or Core
  END_DOC
  character*(32) :: buffer
  buffer = types(t_Simple)
  jast_type = t_Core
  call get_jastrow_jast_type(buffer)
  if (buffer == types(t_Simple)) then
    jast_type = t_Simple
  else if (buffer == types(t_None)) then
    jast_type = t_None
  else if (buffer == types(t_Core)) then
    jast_type = t_Core
  else if (buffer == types(t_Mu)) then
    jast_type = t_Mu
  else
    call abrt(irp_here,'Jastrow type should be (None|Simple|Core|Mu)')
  endif
  call cinfo(irp_here,'jast_type',buffer)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_a_up_up ]
  implicit none
  BEGIN_DOC
! a_{up up} parameters of the Jastrow
  END_DOC
  include '../types.F'
  jast_a_up_up = 0.5
  call get_jastrow_jast_a_up_up(jast_a_up_up)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_a_up_dn ]
  implicit none
  BEGIN_DOC
! a_{up dn} parameters of the Jastrow
  END_DOC
  include '../types.F'
  jast_a_up_dn = 0.5
  call get_jastrow_jast_a_up_dn(jast_a_up_dn)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_b_up_up ]
  implicit none
  BEGIN_DOC
! b_{up up} parameters of the Jastrow
  END_DOC
  include '../types.F'
  jast_b_up_up = 1.
  call get_jastrow_jast_b_up_up(jast_b_up_up)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_b_up_dn ]
  implicit none
  BEGIN_DOC
! b_{up dn} parameters of the Jastrow
  END_DOC
  include '../types.F'
  jast_b_up_dn = 1.
  call get_jastrow_jast_b_up_dn(jast_b_up_dn)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_pen, (nucl_num) ]
  implicit none
  BEGIN_DOC
! penetration parameters of the Jastrow
  END_DOC
  include '../types.F'
  jast_pen(:) = 0.5
  call get_jastrow_jast_pen(jast_pen)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_eeN_e_a, (nucl_num) ]
  implicit none
  BEGIN_DOC
! a parameters of the electron-electron-Nucleus component of the Jastrow
  END_DOC
  include '../types.F'
  jast_eeN_e_a(:) = 0.5
  call get_jastrow_jast_eeN_e_a(jast_eeN_e_a)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_eeN_e_b, (nucl_num) ]
  implicit none
  BEGIN_DOC
! b parameters of the electron-electron-Nucleus component of the Jastrow
  END_DOC
  include '../types.F'
  jast_eeN_e_b(:) = 1.
  call get_jastrow_jast_eeN_e_b(jast_eeN_e_b)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_eeN_N, (nucl_num) ]
  implicit none
  BEGIN_DOC
! penetration parameters of the electron-electron-nucleus component of the Jastrow
  END_DOC
 include '../types.F'
 integer :: i
 jast_eeN_N(:) = 100.
 call get_jastrow_jast_eeN_N(jast_eeN_N)

END_PROVIDER


BEGIN_PROVIDER [ real, jast_core_a1, (nucl_num) ]
  implicit none
  BEGIN_DOC
! parameters of the core Jastrow
  END_DOC
  include '../types.F'
  integer :: i
  do i=1,nucl_num
    if (nucl_charge(i) > 0.) then
      jast_core_a1(i) = 0.6/nucl_charge(i)
    else
      jast_core_a1(i) = 0.
    endif
  enddo
  call get_jastrow_jast_core_a1(jast_core_a1)

END_PROVIDER


 BEGIN_PROVIDER [ real, jast_core_b1, (nucl_num) ]
  implicit none
  BEGIN_DOC
! parameters of the core Jastrow
  END_DOC
  include '../types.F'
  jast_core_b1(:) = max(1.,1. - 0.3 * nucl_charge(:))
  call get_jastrow_jast_core_b1(jast_core_b1)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_mu_erf ]
  implicit none
  include '../types.F'
  jast_mu_erf = 0.5
  call get_jastrow_mu_erf(jast_mu_erf)
END_PROVIDER


BEGIN_PROVIDER  [ integer, jast_1b_type ]
 implicit none
 include '../types.F'
! jast_1b_type = 0 ! no  1body Jastrow
 !jast_1b_type = 2 ! add 1body-Tanh   Jastrow
 !jast_1b_type = 3 ! add 1body-Simple Jastrow
 !jast_1b_type = 4 ! add 1body-RSDFT  Jastrow
 !jast_1b_type = 5 ! add 1body-erf    Jastrow
 jast_1b_type = 6 ! add 1body-Gauss  Jastrow
! call get_jastrow_jast_1b_type(jast_1b_type)
END_PROVIDER

! useful if  jast_1b_type = 2
BEGIN_PROVIDER [ real, jast_1btanh_pen, (nucl_num) ]
  implicit none
  include '../types.F'
  jast_1btanh_pen(:) = 1.0
  call get_jastrow_jast_1btanh_pen(jast_1btanh_pen)
END_PROVIDER

! useful if  jast_1b_type = 5
BEGIN_PROVIDER [ real, jast_1berf_pen, (nucl_num) ]
  implicit none
  include '../types.F'
  jast_1berf_pen(:) = 1.0
  call get_jastrow_jast_1berf_pen(jast_1berf_pen)
END_PROVIDER

! useful if  jast_1b_type = 6
BEGIN_PROVIDER [ real, jast_1bGauss_pen, (nucl_num) ]
  implicit none
  include '../types.F'
  jast_1bGauss_pen(:) = 1.0
  call get_jastrow_jast_1bgauss_pen(jast_1bGauss_pen)
END_PROVIDER
