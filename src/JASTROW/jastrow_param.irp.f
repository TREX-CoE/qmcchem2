! Input data
! ----------

BEGIN_PROVIDER  [ logical, do_jast ]

  BEGIN_DOC  
  ! If true, compute the Jastrow factor
  END_DOC

  implicit none
  include '../types.F'
  do_jast = jast_type /= t_None
  call linfo(irp_here, 'do_jast', do_jast)

END_PROVIDER

! ---

BEGIN_PROVIDER  [ logical, do_jpsi ]

  BEGIN_DOC  
  ! If true, add Jastrow factor to WF
  END_DOC

  implicit none
  include '../types.F'
  do_jpsi = jpsi_type /= t_None
  call linfo(irp_here, 'do_jpsi', do_jpsi)

END_PROVIDER

! ---

BEGIN_PROVIDER [ integer, jast_type ]

  BEGIN_DOC
  ! Type of Jastrow factor : Simple or Core
  END_DOC

  include '../types.F'

  implicit none
  character*(32) :: buffer

  buffer = types(t_None)
  jast_type = t_None
  call get_jastrow_jast_type(buffer)
  if (buffer == types(t_Simple)) then
    jast_type = t_Simple
  else if (buffer == types(t_None)) then
    jast_type = t_None
  else if (buffer == types(t_Core)) then
    jast_type = t_Core
  else if (buffer == types(t_Mu)) then
    jast_type = t_Mu
  else if (buffer == types(t_Mu_1b)) then
    jast_type = t_Mu_1b
  else if (buffer == types(t_Muenv)) then
    jast_type = t_Muenv
  else if (buffer == types(t_Qmckl)) then
    jast_type = t_Qmckl
  else if (buffer == types(t_Mur)) then
    jast_type = t_Mur
    print*, ' do not forget to increase the block time'
  else
    call abrt(irp_here, 'Jastrow type should be (None|Simple|Qmckl|Core|Mu|Mu_1b|Muenv|Mur)')
  endif
  call cinfo(irp_here, 'jast_type',buffer)

END_PROVIDER

! ---

BEGIN_PROVIDER [ integer, jpsi_type ]

  BEGIN_DOC
  ! Type of Jastrow factor used as function not to sample
  END_DOC

  include '../types.F'

  implicit none
  character*(32) :: buffer

  buffer = types(t_Simple)
  jpsi_type = t_Core
  call get_jastrow_jpsi_type(buffer)
  if (buffer == types(t_Simple)) then
    jpsi_type = t_Simple
  else if (buffer == types(t_None)) then
    jpsi_type = t_None
  else if (buffer == types(t_Core)) then
    jpsi_type = t_Core
  else if (buffer == types(t_Mu)) then
    jpsi_type = t_Mu
  else if (buffer == types(t_Mu_1b)) then
    jpsi_type = t_Mu_1b
  else if (buffer == types(t_Muenv)) then
    jpsi_type = t_Muenv
  else if (buffer == types(t_Qmckl)) then
    jpsi_type = t_Qmckl
  else if (buffer == types(t_Mur)) then
    jpsi_type = t_Mur
    print*, ' do not forget to increase the block time'
  else
    call abrt(irp_here, 'jpsi type should be (None|Simple|Qmckl|Core|Mu|Mu_1b|Muenv|Mur)')
  endif
  call cinfo(irp_here, 'jpsi_type',buffer)

END_PROVIDER

! ---

BEGIN_PROVIDER [real, jast_a_up_up]

  BEGIN_DOC
  ! a_{up up} parameters of the Jastrow
  END_DOC

  implicit none

  jast_a_up_up = 0.5
  call get_jastrow_jast_a_up_up(jast_a_up_up)

END_PROVIDER

BEGIN_PROVIDER [real, jast_a_up_dn]

  BEGIN_DOC
  ! a_{up dn} parameters of the Jastrow
  END_DOC

  implicit none

  jast_a_up_dn = 0.5
  call get_jastrow_jast_a_up_dn(jast_a_up_dn)

END_PROVIDER

BEGIN_PROVIDER [real, jast_b_up_up]

  BEGIN_DOC
  ! b_{up up} parameters of the Jastrow
  END_DOC

  implicit none

  jast_b_up_up = 1.
  call get_jastrow_jast_b_up_up(jast_b_up_up)

END_PROVIDER

BEGIN_PROVIDER [real, jast_b_up_dn]

  BEGIN_DOC
  ! b_{up dn} parameters of the Jastrow
  END_DOC

  implicit none

  jast_b_up_dn = 1.
  call get_jastrow_jast_b_up_dn(jast_b_up_dn)

END_PROVIDER

BEGIN_PROVIDER [real, mu_erf]

  BEGIN_DOC
  ! mu parameter
  END_DOC

  implicit none

  mu_erf = 0.5
  call get_ao_two_e_erf_ints_mu_erf(mu_erf)

END_PROVIDER

BEGIN_PROVIDER [real, mu_r_ct]

  BEGIN_DOC
  ! mu_r_ct parameter
  END_DOC

  implicit none

  mu_r_ct = 0.6203504908994001
  call get_tc_keywords_mu_r_ct(mu_r_ct)

END_PROVIDER

! ---

BEGIN_PROVIDER  [ integer, j1b_type ]
  implicit none
  include '../types.F'
  j1b_type = 0 ! no  1body Jastrow
  call get_tc_keywords_j1b_type(j1b_type)
END_PROVIDER

BEGIN_PROVIDER [ real, j1b_pen, (nucl_num) ]
  implicit none
  include '../types.F'
  j1b_pen(:) = 1.0
  call get_tc_keywords_j1b_pen(j1b_pen)
END_PROVIDER

BEGIN_PROVIDER [ real, j1b_pen_coef, (nucl_num) ]
  implicit none
  include '../types.F'
  j1b_pen_coef(:) = 1.0
  call get_tc_keywords_j1b_pen_coef(j1b_pen_coef)
END_PROVIDER

BEGIN_PROVIDER [ double precision, j1b_coeff, (nucl_num) ]
  implicit none
  include '../types.F'
  j1b_coeff(:) = 0.0
  call get_tc_keywords_j1b_coeff(j1b_coeff)
END_PROVIDER

! ---

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
  jast_eeN_e_b(:) = 1.0
  call get_jastrow_jast_eeN_e_b(jast_eeN_e_b)

END_PROVIDER

BEGIN_PROVIDER [ real, jast_eeN_N, (nucl_num) ]
  implicit none
  BEGIN_DOC
! penetration parameters of the electron-electron-nucleus component of the Jastrow
  END_DOC
 include '../types.F'
 integer :: i
 jast_eeN_N(:) = 100.0
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
      jast_core_a1(i) = 0.0
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

