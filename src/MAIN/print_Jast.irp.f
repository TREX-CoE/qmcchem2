program print_Jast

  include '../types.F'
  PROVIDE ezfio_filename

  implicit none
  integer                    :: i, j, k
  double precision           :: accu_tmp, accu_tot, norm
  double precision, external :: qmc_ranf

  print *,  'Number of determinants                   : ', det_num
  print *,  'Number of unique alpha/beta determinants : ', det_alpha_num, det_beta_num
  print *,  'Closed-shell MOs                         : ', mo_closed_num
  print *,  'Number of MOs in determinants            : ', num_present_mos

  print *, ' jastrow mu     ', mu_erf
  print *, ' jastrow 1b type', j1b_type
  print *, ' jastrow 1b     ', j1b_pen, j1b_pen_coef, j1b_coeff
  print *, ' sign of jastrow', sgn_jast

  print *, ' jastrow value:'
  print *, jast_value, jast_value_inv
  print *, ' '

  ! ---

!  accu_tot = 0.d0
!  norm     = 0.d0
!  do i = 1, 1000
!    do k = 1, 3
!      do j = 1, elec_num
!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!      enddo
!    enddo
!    TOUCH elec_coord
!
!    !PROVIDE E_loc
!    !print *,  'E_loc : ', E_loc
!
!    accu_tmp = dabs(jast_Mu_value - jast_Muenv_value)
!
!    accu_tot += accu_tmp
!    norm     += dabs(jast_Muenv_value)
!  enddo
!  print *, ' accu_tot =', accu_tot
!  print *, ' norm     =', norm
!  print *, ' '
!
!  ! ---
!
!!  accu_tot = 0.d0
!!  norm     = 0.d0
!!  do i = 1, 1000
!!    do k = 1, 3
!!      do j = 1, elec_num
!!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!!      enddo
!!    enddo
!!    TOUCH elec_coord
!!
!!    accu_tmp = 0.d0
!!    do j = 1, elec_num
!!      accu_tmp += dabs(jast_grad_jast_inv_x(j) - jast_elec_Muenv_grad_x(j))
!!      accu_tmp += dabs(jast_grad_jast_inv_y(j) - jast_elec_Muenv_grad_y(j))
!!      accu_tmp += dabs(jast_grad_jast_inv_z(j) - jast_elec_Muenv_grad_z(j))
!!      norm     += dabs(jast_elec_Muenv_grad_x(j))
!!      norm     += dabs(jast_elec_Muenv_grad_y(j))
!!      norm     += dabs(jast_elec_Muenv_grad_z(j))
!!    enddo
!!
!!    accu_tot += accu_tmp
!!  enddo
!!  print *, ' accu_tot grad =', accu_tot
!!  print *, ' norm          =', norm
!
!  accu_tot = 0.d0
!  norm     = 0.d0
!  do i = 1, 1000
!
!    do k = 1, 3
!      do j = 1, elec_num
!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!      enddo
!    enddo
!    TOUCH elec_coord
!
!    accu_tmp = 0.d0
!    do j = 1, elec_num
!      accu_tmp += dabs(jast_elec_Mu_grad_x(j) - jast_elec_Muenv_grad_x(j))
!      accu_tmp += dabs(jast_elec_Mu_grad_y(j) - jast_elec_Muenv_grad_y(j))
!      accu_tmp += dabs(jast_elec_Mu_grad_z(j) - jast_elec_Muenv_grad_z(j))
!      norm     += dabs(jast_elec_Muenv_grad_x(j))
!      norm     += dabs(jast_elec_Muenv_grad_y(j))
!      norm     += dabs(jast_elec_Muenv_grad_z(j))
!    enddo
!
!    accu_tot += accu_tmp
!  enddo
!  print *, ' accu_tot grad =', accu_tot
!  print *, ' norm          =', norm
!  print *, ' '
!
!  ! ---
!
!!  accu_tot = 0.d0
!!  norm     = 0.d0
!!  do i = 1, 1000
!!    do k = 1, 3
!!      do j = 1, elec_num
!!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!!      enddo
!!    enddo
!!    TOUCH elec_coord
!!
!!    accu_tmp = 0.d0
!!    do j = 1, elec_num
!!      accu_tmp += dabs(jast_lapl1(j) - jast_elec_Muenv_lapl(j))
!!      norm     += dabs(jast_elec_Muenv_lapl(j))
!!    enddo
!!
!!    accu_tot += accu_tmp
!!  enddo
!!  print *, ' accu_tot lapl =', accu_tot
!!  print *, ' norm          =', norm
!
!  accu_tot = 0.d0
!  norm     = 0.d0
!  do i = 1, 1000
!
!    do k = 1, 3
!      do j = 1, elec_num
!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!      enddo
!    enddo
!    TOUCH elec_coord
!
!    accu_tmp = 0.d0
!    do j = 1, elec_num
!      accu_tmp += dabs(jast_elec_Mu_lapl(j) - jast_elec_Muenv_lapl(j))
!      norm     += dabs(jast_elec_Muenv_lapl(j))
!    enddo
!
!    accu_tot += accu_tmp
!  enddo
!  print *, ' accu_tot lapl =', accu_tot
!  print *, ' norm          =', norm
!
!  ! ---


  print*, ' jast_type', jast_type
  if(jast_type .eq. t_Muenv) then
    call check_Muenv_grad()
    call check_Muenv_lapl()
  elseif(jast_type .eq. t_Mur) then
    !call debug_hf_density()
    call check_Mur_grad()
    call check_Mur_lapl()
  endif

end

! ---

BEGIN_TEMPLATE
subroutine check_$X_grad()

  implicit none
  integer                    :: i, j, k
  double precision           :: accu_tmp, accu_tot, norm, thr
  double precision, external :: qmc_ranf

  thr = 1d-6

  accu_tot = 0.d0
  norm     = 0.d0
!  do i = 1, 1000
!
    do k = 1, 3
      do j = 1, elec_num
        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
      enddo
    enddo
    TOUCH elec_coord

    do j = 1, elec_num

      accu_tmp = dabs(jast_elec_$X_grad_x(j) - jast_elec_$X_grad_x_num(j))
      if(accu_tmp .gt. thr) then
        print *, ' problem on ', j
        print *, ' jast_elec_$X_grad_x    ', jast_elec_$X_grad_x    (j)
        print *, ' jast_elec_$X_grad_x_num', jast_elec_$X_grad_x_num(j)
      endif
      accu_tot += accu_tmp
      norm     += dabs(jast_elec_$X_grad_x(j))

      accu_tmp = dabs(jast_elec_$X_grad_y(j) - jast_elec_$X_grad_y_num(j))
      if(accu_tmp .gt. thr) then
        print *, ' problem on ', j
        print *, ' jast_elec_$X_grad_y    ', jast_elec_$X_grad_y    (j)
        print *, ' jast_elec_$X_grad_y_num', jast_elec_$X_grad_y_num(j)
      endif
      accu_tot += accu_tmp
      norm     += dabs(jast_elec_$X_grad_y(j))

      accu_tmp = dabs(jast_elec_$X_grad_z(j) - jast_elec_$X_grad_z_num(j))
      if(accu_tmp .gt. thr) then
        print *, ' problem on ', j
        print *, ' jast_elec_$X_grad_z    ', jast_elec_$X_grad_z    (j)
        print *, ' jast_elec_$X_grad_z_num', jast_elec_$X_grad_z_num(j)
      endif
      accu_tot += accu_tmp
      norm     += dabs(jast_elec_$X_grad_z(j))

    enddo
!  enddo

  print *, ' accu_tot =', accu_tot
  print *, ' norm     =', norm
  print *, ' '

end subroutine check_$X_grad

subroutine check_$X_lapl()

  implicit none
  integer                    :: i, j, k
  double precision           :: accu_tmp, accu_tot, norm, thr
  double precision, external :: qmc_ranf

  thr = 1d-6

  accu_tot = 0.d0
  norm     = 0.d0
!  do i = 1, 1000
!
    do k = 1, 3
      do j = 1, elec_num
        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
      enddo
    enddo
    TOUCH elec_coord

    do j = 1, elec_num

      accu_tmp = dabs(jast_elec_$X_lapl(j) - jast_elec_$X_lapl_num(j))
      if(accu_tmp .gt. thr) then
        print *, ' problem on ', j
        print *, ' jast_elec_$X_lapl    ', jast_elec_$X_lapl    (j)
        print *, ' jast_elec_$X_lapl_num', jast_elec_$X_lapl_num(j)
      endif
      accu_tot += accu_tmp
      norm     += dabs(jast_elec_$X_lapl(j))

    enddo
!  enddo

  print *, ' accu_tot =', accu_tot
  print *, ' norm     =', norm
  print *, ' '

end subroutine check_$X_lapl
SUBST [X]
Muenv ;;
Mur   ;;
END_TEMPLATE


subroutine debug_hf_density()

  implicit none
  integer                    :: mu, i, j, k
  real                       :: r(3)
  real                       :: ao_val1, ao_der1(3), ao_lap1
  real                       :: ao_val2, ao_der2(3), ao_lap2
  real                       :: accu_tmp, accu_tot, norm, thr
  double precision, external :: qmc_ranf

  print*, ' ** debug HF density **'

  thr = 1e-5

  accu_tot = 0.0
  norm     = 0.0
  do i = 1, 10

    r(1) = 1.5 * (0.5-qmc_ranf())
    r(2) = 1.5 * (0.5-qmc_ranf())
    r(3) = 1.5 * (0.5-qmc_ranf())

    do mu = 1, ao_num

      call get_ao_val_der_lap(mu, r, ao_val1, ao_der1, ao_lap1)

      point(1) = r(1) 
      point(2) = r(2) 
      point(3) = r(3) 
      TOUCH point
      ao_val2    = ao_value_p(mu)
      ao_der2(1) = ao_grad_p (mu,1)
      ao_der2(2) = ao_grad_p (mu,2)
      ao_der2(3) = ao_grad_p (mu,3)
      ao_lap2    = ao_lapl_p (mu)

      accu_tmp = abs(ao_val2 - ao_val1)
      if(accu_tmp .gt. thr) then
        print *, ' problem on val on', mu
        print *, ao_val1, ao_val2
      endif
      accu_tot += accu_tmp
      norm     += abs(ao_val2)

      accu_tmp = abs(ao_der2(1) - ao_der1(1))
      if(accu_tmp .gt. thr) then
        print *, ' problem on der_x on', mu
        print *, ao_der1(1), ao_der2(1)
      endif
      accu_tot += accu_tmp
      norm     += abs(ao_der2(1))

      accu_tmp = abs(ao_der2(2) - ao_der1(2))
      if(accu_tmp .gt. thr) then
        print *, ' problem on der_y on', mu
        print *, ao_der1(2), ao_der2(2)
      endif
      accu_tot += accu_tmp
      norm     += abs(ao_der2(2))

      accu_tmp = abs(ao_der2(3) - ao_der1(3))
      if(accu_tmp .gt. thr) then
        print *, ' problem on der_z on', mu
        print *, ao_der1(3), ao_der2(3)
      endif
      accu_tot += accu_tmp
      norm     += abs(ao_der2(3))

      accu_tmp = abs(ao_lap2 - ao_lap1)
      if(accu_tmp .gt. thr) then
        print *, ' problem on lap on', mu
        print *, ao_lap1, ao_lap2
      endif
      accu_tot += accu_tmp
      norm     += abs(ao_lap2)

    enddo
  enddo

  print *, ' accu_tot =', accu_tot
  print *, ' norm     =', norm

  return
end subroutine debug_hf_density


