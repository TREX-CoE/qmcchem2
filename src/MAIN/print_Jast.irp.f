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
  print *, ' jastrow 1b     ', j1b_pen, j1b_coeff
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
!    accu_tmp = dabs(jast_Mu_value - jast_Mu_env3_value)
!
!    accu_tot += accu_tmp
!    norm     += dabs(jast_Mu_env3_value)
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
!!      accu_tmp += dabs(jast_grad_jast_inv_x(j) - jast_elec_Mu_env3_grad_x(j))
!!      accu_tmp += dabs(jast_grad_jast_inv_y(j) - jast_elec_Mu_env3_grad_y(j))
!!      accu_tmp += dabs(jast_grad_jast_inv_z(j) - jast_elec_Mu_env3_grad_z(j))
!!      norm     += dabs(jast_elec_Mu_env3_grad_x(j))
!!      norm     += dabs(jast_elec_Mu_env3_grad_y(j))
!!      norm     += dabs(jast_elec_Mu_env3_grad_z(j))
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
!      accu_tmp += dabs(jast_elec_Mu_grad_x(j) - jast_elec_Mu_env3_grad_x(j))
!      accu_tmp += dabs(jast_elec_Mu_grad_y(j) - jast_elec_Mu_env3_grad_y(j))
!      accu_tmp += dabs(jast_elec_Mu_grad_z(j) - jast_elec_Mu_env3_grad_z(j))
!      norm     += dabs(jast_elec_Mu_env3_grad_x(j))
!      norm     += dabs(jast_elec_Mu_env3_grad_y(j))
!      norm     += dabs(jast_elec_Mu_env3_grad_z(j))
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
!!      accu_tmp += dabs(jast_lapl1(j) - jast_elec_Mu_env3_lapl(j))
!!      norm     += dabs(jast_elec_Mu_env3_lapl(j))
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
!      accu_tmp += dabs(jast_elec_Mu_lapl(j) - jast_elec_Mu_env3_lapl(j))
!      norm     += dabs(jast_elec_Mu_env3_lapl(j))
!    enddo
!
!    accu_tot += accu_tmp
!  enddo
!  print *, ' accu_tot lapl =', accu_tot
!  print *, ' norm          =', norm
!
!  ! ---


  PROVIDE jast_type
  print*, ' jast_type', jast_type
  if(jast_type .eq. t_Mu_env2) then
    call check_Mu_env2_grad()
    call check_Mu_env2_lapl()
  elseif(jast_type .eq. t_Mu_env3) then
    call check_Mu_env3_grad()
    call check_Mu_env3_lapl()
  elseif(jast_type .eq. t_Mu_env5) then
    call check_Mu_env5_grad()
    call check_Mu_env5_lapl()
  endif

end



! ---

!subroutine check_jmu_env3_grad()
!
!  implicit none
!  integer                    :: i, j, k
!  double precision           :: accu_tmp, accu_tot, norm, thr
!  double precision, external :: qmc_ranf
!
!  thr = 1d-6
!
!  accu_tot = 0.d0
!  norm     = 0.d0
!!  do i = 1, 1000
!!
!    do k = 1, 3
!      do j = 1, elec_num
!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!      enddo
!    enddo
!    TOUCH elec_coord
!
!    do j = 1, elec_num
!
!      accu_tmp = dabs(jast_elec_Mu_env3_grad_x(j) - jast_elec_Mu_env3_grad_x_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env3_grad_x    ', jast_elec_Mu_env3_grad_x    (j)
!        print *, ' jast_elec_Mu_env3_grad_x_num', jast_elec_Mu_env3_grad_x_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env3_grad_x(j))
!
!      accu_tmp = dabs(jast_elec_Mu_env3_grad_y(j) - jast_elec_Mu_env3_grad_y_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env3_grad_y    ', jast_elec_Mu_env3_grad_y    (j)
!        print *, ' jast_elec_Mu_env3_grad_y_num', jast_elec_Mu_env3_grad_y_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env3_grad_y(j))
!
!      accu_tmp = dabs(jast_elec_Mu_env3_grad_z(j) - jast_elec_Mu_env3_grad_z_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env3_grad_z    ', jast_elec_Mu_env3_grad_z    (j)
!        print *, ' jast_elec_Mu_env3_grad_z_num', jast_elec_Mu_env3_grad_z_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env3_grad_z(j))
!
!    enddo
!!  enddo
!
!  print *, ' accu_tot =', accu_tot
!  print *, ' norm     =', norm
!  print *, ' '
!
!end subroutine check_jmu_env3_grad
!
!! ---
!
!subroutine check_jmu_env3_lapl()
!
!  implicit none
!  integer                    :: i, j, k
!  double precision           :: accu_tmp, accu_tot, norm, thr
!  double precision, external :: qmc_ranf
!
!  thr = 1d-6
!
!  accu_tot = 0.d0
!  norm     = 0.d0
!!  do i = 1, 1000
!!
!    do k = 1, 3
!      do j = 1, elec_num
!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!      enddo
!    enddo
!    TOUCH elec_coord
!
!    do j = 1, elec_num
!
!      accu_tmp = dabs(jast_elec_Mu_env3_lapl(j) - jast_elec_Mu_env3_lapl_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env3_lapl    ', jast_elec_Mu_env3_lapl    (j)
!        print *, ' jast_elec_Mu_env3_lapl_num', jast_elec_Mu_env3_lapl_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env3_lapl(j))
!
!    enddo
!!  enddo
!
!  print *, ' accu_tot =', accu_tot
!  print *, ' norm     =', norm
!  print *, ' '
!
!end subroutine check_jmu_env3_lapl
!
!! ---
!
!subroutine check_jmu_env5_grad()
!
!  implicit none
!  integer                    :: i, j, k
!  double precision           :: accu_tmp, accu_tot, norm, thr
!  double precision, external :: qmc_ranf
!
!  thr = 1d-6
!
!  accu_tot = 0.d0
!  norm     = 0.d0
!!  do i = 1, 1000
!!
!    do k = 1, 3
!      do j = 1, elec_num
!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!      enddo
!    enddo
!    TOUCH elec_coord
!
!    do j = 1, elec_num
!
!      accu_tmp = dabs(jast_elec_Mu_env5_grad_x(j) - jast_elec_Mu_env5_grad_x_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env5_grad_x    ', jast_elec_Mu_env5_grad_x    (j)
!        print *, ' jast_elec_Mu_env5_grad_x_num', jast_elec_Mu_env5_grad_x_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env5_grad_x(j))
!
!      accu_tmp = dabs(jast_elec_Mu_env5_grad_y(j) - jast_elec_Mu_env5_grad_y_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env5_grad_y    ', jast_elec_Mu_env5_grad_y    (j)
!        print *, ' jast_elec_Mu_env5_grad_y_num', jast_elec_Mu_env5_grad_y_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env5_grad_y(j))
!
!      accu_tmp = dabs(jast_elec_Mu_env5_grad_z(j) - jast_elec_Mu_env5_grad_z_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env5_grad_z    ', jast_elec_Mu_env5_grad_z    (j)
!        print *, ' jast_elec_Mu_env5_grad_z_num', jast_elec_Mu_env5_grad_z_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env5_grad_z(j))
!
!    enddo
!!  enddo
!
!  print *, ' accu_tot =', accu_tot
!  print *, ' norm     =', norm
!  print *, ' '
!
!end subroutine check_jmu_env5_grad
!
!! ---
!
!subroutine check_jmu_env5_lapl()
!
!  implicit none
!  integer                    :: i, j, k
!  double precision           :: accu_tmp, accu_tot, norm, thr
!  double precision, external :: qmc_ranf
!
!  thr = 1d-6
!
!  accu_tot = 0.d0
!  norm     = 0.d0
!!  do i = 1, 1000
!!
!    do k = 1, 3
!      do j = 1, elec_num
!        elec_coord(j,k) += 1.5 * (0.5-qmc_ranf())
!      enddo
!    enddo
!    TOUCH elec_coord
!
!    do j = 1, elec_num
!
!      accu_tmp = dabs(jast_elec_Mu_env5_lapl(j) - jast_elec_Mu_env5_lapl_num(j))
!      if(accu_tmp .gt. thr) then
!        print *, ' problem on ', j
!        print *, ' jast_elec_Mu_env5_lapl    ', jast_elec_Mu_env5_lapl    (j)
!        print *, ' jast_elec_Mu_env5_lapl_num', jast_elec_Mu_env5_lapl_num(j)
!      endif
!      accu_tot += accu_tmp
!      norm     += dabs(jast_elec_Mu_env5_lapl(j))
!
!    enddo
!!  enddo
!
!  print *, ' accu_tot =', accu_tot
!  print *, ' norm     =', norm
!  print *, ' '
!
!end subroutine check_jmu_env5_lapl

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
Mu_env2 ;;
Mu_env3 ;;
Mu_env5 ;;
END_TEMPLATE

