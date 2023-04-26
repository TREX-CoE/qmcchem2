program print_he

  PROVIDE ezfio_filename

  implicit none
  integer                    :: i, n_theta
  double precision           :: pi, d_theta, theta 

  print *,  'Number of determinants                   : ', det_num
  print *,  'Number of unique alpha/beta determinants : ', det_alpha_num, det_beta_num
  print *,  'Closed-shell MOs                         : ', mo_closed_num
  print *,  'Number of MOs in determinants            : ', num_present_mos

  print *, ' jastrow mu     ', mu_erf
  print *, ' jastrow 1b type', j1b_type
  print *, ' jastrow 1b     ', j1b_pen, j1b_coeff
  print *, ' sign of jastrow', sgn_jast

  print *, ' jastrow value:'
  print *, jast_value, jast_value_inv
  print *, ' '

  ! ---

  pi      = 3.14d0
  n_theta = 250
  d_theta = 2.d0 * pi / dble(n_theta)

  do i = 1, n_theta

    theta = -pi + dble(i-1) * d_theta

    elec_coord(1,1) = 0.5d0
    elec_coord(1,2) = 0.0d0
    elec_coord(1,3) = 0.0d0
    elec_coord(2,1) = 0.5d0 * dcos(theta)
    elec_coord(2,2) = 0.5d0 * dsin(theta)
    elec_coord(2,3) = 0.0d0
    TOUCH elec_coord

    print *, theta, psidet_right_value, psi_value

  enddo


end

