
! ---

BEGIN_PROVIDER [ real, mo_left_coef_input, (ao_num_8,mo_tot_num) ]

  BEGIN_DOC
  ! left molecular orbital coefficients read from the input file
  END_DOC

  implicit none
  integer           :: i, j
  real, allocatable :: buffer(:,:)

  allocate (buffer(ao_num,mo_tot_num))

  buffer = 0.d0
  call get_bi_ortho_mos_mo_l_coef(buffer)
  do i = 1, mo_tot_num
    do j = 1, ao_num
      mo_left_coef_input(j,i) = buffer(j,i)
    enddo
    call set_order(mo_left_coef_input(1,i), ao_nucl_sort_idx, ao_num)
    do j = ao_num+1, ao_num_8
      mo_left_coef_input(j,i) = 0.d0
    enddo
  enddo

  deallocate(buffer)

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, mo_left_coef, (ao_num_8,mo_num_8) ]

  BEGIN_DOC
  !left mMolecular orbital coefficients
  END_DOC

  implicit none
  integer :: i, j

  do j = 1, mo_num
    do i = 1, ao_num_8
      mo_left_coef(i,j) = mo_left_coef_input(i,j)
    enddo
  enddo
  do j = mo_num+1, mo_num_8
    !DIR$ VECTOR ALIGNED
    do i = 1, ao_num_8
      mo_left_coef(i,j) = 0.d0
    enddo
  enddo

  ! Input MOs are not needed any more
  FREE mo_left_coef_input

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, mo_left_coef_transp, (mo_num_8,ao_num_8) ]

  BEGIN_DOC
  ! Transpose of the left molecular orbital coefficients
  END_DOC

  implicit none
  call transpose(mo_left_coef, ao_num_8, mo_left_coef_transp, mo_num_8, ao_num_8, mo_num_8)

END_PROVIDER

! ---

 BEGIN_PROVIDER [ integer, mo_left_coef_transp_non_zero_idx, (0:mo_num,ao_num) ]
&BEGIN_PROVIDER [ real, mo_left_coef_transp_sparsity ]

  BEGIN_DOC
  ! Indices of the non-zero elements of the transpose of the left molecular
  ! orbital coefficients
  END_DOC

  implicit none
  integer :: i, j
  integer :: idx

  mo_left_coef_transp_sparsity = 0.
  do j = 1, ao_num
    idx = 0
    do i = 1, mo_num
      if(mo_left_coef_transp(i,j) /= 0.) then
        idx += 1
        mo_left_coef_transp_non_zero_idx(idx,j) = i
      endif
    enddo
    mo_left_coef_transp_non_zero_idx(0,j) = idx
    mo_left_coef_transp_sparsity += float(idx)
  enddo

  mo_left_coef_transp_sparsity *= 1./(mo_num*ao_num)

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, mo_left_coef_transp_present, (num_present_mos_8,ao_num_8) ]

  BEGIN_DOC
  ! mo_left_coef_transp without MOs absent in all determinants
  END_DOC

  implicit none
  integer :: i,j,n

  mo_left_coef_transp_present = 0.
  do i = 1, ao_num
    do j = 1, num_present_mos
      mo_left_coef_transp_present(j,i) = mo_left_coef_transp(present_mos(j),i)
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ real, mo_left_value_transp,  (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_left_grad_transp_x, (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_left_grad_transp_y, (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_left_grad_transp_z, (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_left_lapl_transp, (mo_num_8,elec_num) ]

  BEGIN_DOC
  ! Values, gradients, laplacians of the left molecular orbitals
  ! Arrays are padded for efficiency
  END_DOC

  implicit none
  integer :: i, j, k, l, m

  PROVIDE primitives_reduced

  if (do_nucl_fitcusp) then
    PROVIDE nucl_left_fitcusp_param
    PROVIDE nucl_elec_dist_vec
    PROVIDE nucl_elec_dist_inv
    mo_left_value_transp  = 0.
    mo_left_grad_transp_x = 0.
    mo_left_grad_transp_y = 0.
    mo_left_grad_transp_z = 0.
    mo_left_lapl_transp   = 0.
  endif

  if (use_qmckl) then

    print*, ' NOT IMPLEMENTED YET'
    stop

  else

    do i = 1, elec_num

      if(i>1) then
        ao_elec = i
        TOUCH ao_elec
      endif

      if(num_present_mos == mo_num) then
        call sparse_full_mv(mo_left_coef_transp,mo_num_8,              &
            ao_value_block(1),ao_num_8,                                &
            ao_grad_block_x(1),                                        &
            ao_grad_block_y(1),                                        &
            ao_grad_block_z(1),                                        &
            ao_lapl_block(1),                                          &
            ao_value_non_zero_idx(0),                                  &
            mo_left_value_transp(1,i),mo_num_8,                        &
            mo_left_grad_transp_x(1,i),                                &
            mo_left_grad_transp_y(1,i),                                &
            mo_left_grad_transp_z(1,i),                                &
            mo_left_lapl_transp(1,i),                                  &
            ao_num)

      else
        call sparse_full_mv(mo_left_coef_transp_present,num_present_mos_8,  &
            ao_value_block(1),ao_num_8,                                     &
            ao_grad_block_x(1),                                             &
            ao_grad_block_y(1),                                             &
            ao_grad_block_z(1),                                             &
            ao_lapl_block(1),                                               &
            ao_value_non_zero_idx(0),                                       &
            mo_left_value_transp(1,i),mo_num_8,                             &
            mo_left_grad_transp_x(1,i),                                     &
            mo_left_grad_transp_y(1,i),                                     &
            mo_left_grad_transp_z(1,i),                                     &
            mo_left_lapl_transp(1,i),                                       &
            ao_num)

      endif
    enddo ! i

    ao_elec = 1
    SOFT_TOUCH ao_elec

  endif

  if(do_nucl_fitcusp) then
    real :: r, r2, r_inv, d, expzr, Z, Z2, a, b, c, phi, rx, ry, rz
    do i=1,elec_num
      do k = 1, nucl_num
        r = nucl_elec_dist(k,i)
        if (r > nucl_fitcusp_radius(k)) then
          cycle
        endif
        r_inv = nucl_elec_dist_inv(k,i)
        do j = 1, num_present_mos
          mo_left_value_transp(j,i) = mo_left_value_transp(j,i) + nucl_left_fitcusp_param(1,j,k) +&
              r * (nucl_left_fitcusp_param(2,j,k) +                       &
              r * (nucl_left_fitcusp_param(3,j,k) +                       &
              r *  nucl_left_fitcusp_param(4,j,k) ))
          mo_left_lapl_transp(j,i) = mo_left_lapl_transp(j,i) +                &
              nucl_left_fitcusp_param(2,j,k)*(r_inv+r_inv) +              &
              6.*nucl_left_fitcusp_param(3,j,k) +                         &
              r * 12.*nucl_left_fitcusp_param(4,j,k)
          c = r_inv * (nucl_left_fitcusp_param(2,j,k) +                   &
              r * (2.*nucl_left_fitcusp_param(3,j,k) +                    &
              r *  3.*nucl_left_fitcusp_param(4,j,k) ))
          mo_left_grad_transp_x(j,i) = mo_left_grad_transp_x(j,i) + nucl_elec_dist_vec(1,k,i)*c
          mo_left_grad_transp_y(j,i) = mo_left_grad_transp_y(j,i) + nucl_elec_dist_vec(2,k,i)*c
          mo_left_grad_transp_z(j,i) = mo_left_grad_transp_z(j,i) + nucl_elec_dist_vec(3,k,i)*c
        enddo
        ! It is safe to exit here because a core electron is close to only one nucleus
        exit
      enddo ! k
    enddo ! i

  endif

  if(do_prepare) then
    real                           :: lambda, t
    ! Scale off-diagonal elements
    t = prepare_walkers_t
    do i = 1, mo_num
      do j = 1, elec_alpha_num
        if (i /= j) then
          mo_left_value_transp (i,j) *= t
          mo_left_grad_transp_x(i,j) *= t
          mo_left_grad_transp_y(i,j) *= t
          mo_left_grad_transp_z(i,j) *= t
          mo_left_lapl_transp  (i,j) *= t
        endif
      enddo
      do j = 1, elec_beta_num
        if (i /= j) then
          mo_left_value_transp (i,j+elec_alpha_num) *= t
          mo_left_grad_transp_x(i,j+elec_alpha_num) *= t
          mo_left_grad_transp_y(i,j+elec_alpha_num) *= t
          mo_left_grad_transp_z(i,j+elec_alpha_num) *= t
          mo_left_lapl_transp  (i,j+elec_alpha_num) *= t
        endif
      enddo
    enddo
  endif

  if (num_present_mos < mo_num) then
    do i=1,elec_num
      do j=num_present_mos,1,-1
        mo_left_value_transp (present_mos(j),i) = mo_left_value_transp (j,i)
        mo_left_grad_transp_x(present_mos(j),i) = mo_left_grad_transp_x(j,i)
        mo_left_grad_transp_y(present_mos(j),i) = mo_left_grad_transp_y(j,i)
        mo_left_grad_transp_z(present_mos(j),i) = mo_left_grad_transp_z(j,i)
        mo_left_lapl_transp  (present_mos(j),i) = mo_left_lapl_transp  (j,i)
        if (present_mos(j) == j) then
          exit
        endif
      enddo
    enddo
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, mo_left_value, (elec_num_8,mo_num) ]

  BEGIN_DOC
  ! Values of the molecular orbitals
  END_DOC

  implicit none
  integer       :: i, j
  integer, save :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    PROVIDE primitives_reduced
    !DIR$ VECTOR ALIGNED
    mo_left_value = 0.
  endif

  call transpose(mo_left_value_transp(1,1), mo_num_8, mo_left_value, elec_num_8, mo_num, elec_num)

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, mo_left_grad_x, (elec_num_8,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_left_grad_y, (elec_num_8,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_left_grad_z, (elec_num_8,mo_num) ]

  BEGIN_DOC
  ! Gradients of the molecular orbitals
  END_DOC

  implicit none
  integer       :: i, j
  integer, save :: ifirst = 0

  if (ifirst == 0) then
    !DIR$ VECTOR ALIGNED
    mo_left_grad_x = 0.d0
    !DIR$ VECTOR ALIGNED
    mo_left_grad_y = 0.d0
    !DIR$ VECTOR ALIGNED
    mo_left_grad_z = 0.d0
    ifirst = 1
    PROVIDE primitives_reduced
  endif
  ! Transpose x last for cache efficiency
  call transpose_to_dp(mo_left_grad_transp_y(1,1), mo_num_8, mo_left_grad_y(1,1), elec_num_8, mo_num, elec_num)
  call transpose_to_dp(mo_left_grad_transp_z(1,1), mo_num_8, mo_left_grad_z(1,1), elec_num_8, mo_num, elec_num)
  call transpose_to_dp(mo_left_grad_transp_x(1,1), mo_num_8, mo_left_grad_x(1,1), elec_num_8, mo_num, elec_num)

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, mo_left_lapl, (elec_num_8,mo_num) ]

  BEGIN_DOC
  ! Laplacians of the molecular orbitals
  END_DOC

  implicit none
  integer       :: i, j
  integer, save :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    PROVIDE primitives_reduced
    !DIR$ VECTOR ALIGNED
    mo_left_lapl = 0.d0
  endif
  call transpose_to_dp(mo_left_lapl_transp(1,1), mo_num_8, mo_left_lapl, elec_num_8, mo_num, elec_num)

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, mo_left_grad_lapl_alpha, (4,elec_alpha_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_left_grad_lapl_beta , (4,elec_alpha_num+1:elec_num,mo_num) ]

  BEGIN_DOC
  ! Gradients and laplacian
  END_DOC
 
  implicit none
  integer :: i,j
  do j=1,mo_num
     do i=1,elec_alpha_num
       mo_left_grad_lapl_alpha(1,i,j) = mo_left_grad_transp_x(j,i)
       mo_left_grad_lapl_alpha(2,i,j) = mo_left_grad_transp_y(j,i)
       mo_left_grad_lapl_alpha(3,i,j) = mo_left_grad_transp_z(j,i)
       mo_left_grad_lapl_alpha(4,i,j) = mo_left_lapl_transp  (j,i)
     enddo
   enddo
   do j=1,mo_num
     do i=elec_alpha_num+1,elec_num
       mo_left_grad_lapl_beta(1,i,j) = mo_left_grad_transp_x(j,i)
       mo_left_grad_lapl_beta(2,i,j) = mo_left_grad_transp_y(j,i)
       mo_left_grad_lapl_beta(3,i,j) = mo_left_grad_transp_z(j,i)
       mo_left_grad_lapl_beta(4,i,j) = mo_left_lapl_transp  (j,i)
   enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, mo_left_grad_lapl_transp, (4,mo_num,elec_num) ]

  BEGIN_DOC
  ! Gradients and laplacian
  END_DOC

  implicit none
  integer :: i,j
  do i=1,elec_num
    do j=1,mo_num
      mo_left_grad_lapl_transp(1,j,i) = mo_left_grad_transp_x(j,i)
      mo_left_grad_lapl_transp(2,j,i) = mo_left_grad_transp_y(j,i)
      mo_left_grad_lapl_transp(3,j,i) = mo_left_grad_transp_z(j,i)
      mo_left_grad_lapl_transp(4,j,i) = mo_left_lapl_transp  (j,i)
    enddo
  enddo

END_PROVIDER

! ---

!-----------------
! Fit cusp
!-----------------


BEGIN_PROVIDER [ double precision , mo_left_value_at_nucl, (mo_num_8,nucl_num) ]

  BEGIN_DOC
  ! Values of the molecular orbitals at the nucleus without the
  ! S components of the current nucleus
  END_DOC

  implicit none
  integer :: i, j, k, l
  real    :: ao_value_at_nucl_no_S(ao_num)

  PROVIDE mo_left_fitcusp_normalization_before
  do k = 1, nucl_num
    point(1) = nucl_coord(k,1)
    point(2) = nucl_coord(k,2)
    point(3) = nucl_coord(k,3)
    TOUCH point

    PROVIDE ao_value_p

    do i = 1, ao_num
      if(ao_nucl(i) /= k) then
        ao_value_at_nucl_no_S(i) = ao_value_p(i)
      else
        ao_value_at_nucl_no_S(i) = 0.
      endif
    enddo

    integer :: jj
    do jj = 1, num_present_mos
      j = present_mos(jj)
      mo_left_value_at_nucl(j,k) = 0.
      !DIR$ VECTOR ALIGNED
      do i = 1, ao_num
        mo_left_value_at_nucl(j,k) = mo_left_value_at_nucl(j,k) + mo_left_coef(i,j) * ao_value_at_nucl_no_S(i)
      enddo
    enddo

  enddo
  FREE ao_value_p ao_grad_p ao_lapl_p ao_axis_grad_p ao_oned_grad_p ao_oned_prim_grad_p
  FREE ao_oned_lapl_p ao_axis_lapl_p ao_oned_prim_lapl_p ao_oned_p ao_oned_prim_p
  FREE ao_axis_p ao_axis_power_p
  SOFT_TOUCH point

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, mo_left_fitcusp_normalization_before, (mo_tot_num) ]

  BEGIN_DOC
  ! Renormalization factor of MOs due to cusp fitting
  END_DOC

  include 'constants.F'

  implicit none
  integer          :: i, j, k, l
  double precision :: dr, r, f, t
  integer, save    :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    mo_left_fitcusp_normalization_before = 0.d0
    do k = 1, nucl_num
      dr = nucl_fitcusp_radius(k)*1.d-2
      point(1) = nucl_coord(k,1)
      point(2) = nucl_coord(k,2)
      point(3) = nucl_coord(k,3)-dr
      do l = 1, 101
        r = point(3) + dr
        point(3) = r
        TOUCH point
        f = dfour_pi*r*r*dr
        do i = 1, mo_tot_num
          mo_left_fitcusp_normalization_before(i) += f * mo_left_value_p(i)**2
        enddo
      enddo
    enddo
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, mo_left_fitcusp_normalization_after, (mo_tot_num) ]

  BEGIN_DOC
  ! Renormalization factor of MOs due to cusp fitting
  END_DOC

  include 'constants.F'

  implicit none
  integer          :: i, j, k, l
  double precision :: dr, r, f, t, t2
  integer, save    :: ifirst = 0

  PROVIDE primitives_reduced
  if (ifirst == 0) then
    ifirst = 1
    mo_left_fitcusp_normalization_after = 0.d0
    do k = 1, nucl_num
      dr = nucl_fitcusp_radius(k)*1.d-2
      point(1) = nucl_coord(k,1)
      point(2) = nucl_coord(k,2)
      point(3) = nucl_coord(k,3) - dr
      do l=1,101
        point(3) = point(3)+ dr
        TOUCH point nucl_left_fitcusp_param primitives_reduced mo_left_coef
        r = point(3)
        f = dfour_pi*r*r*dr
        do i = 1, mo_num
          t = 0.d0
          do j = 1, ao_num
            if ( (ao_nucl(j) /= k).or.(ao_power(j,4) > 0) ) then
              t = t + mo_left_coef(j,i) * ao_value_p(j)
            endif
          enddo
          t = t +  nucl_left_fitcusp_param(1,i,k) +                       &
              r * (nucl_left_fitcusp_param(2,i,k) +                       &
              r * (nucl_left_fitcusp_param(3,i,k) +                       &
              r *  nucl_left_fitcusp_param(4,i,k) ))
          mo_left_fitcusp_normalization_after(i) += t*t*f
        enddo
      enddo
    enddo
  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, mo_left_value_at_fitcusp_radius, (mo_num_8,nucl_num) ]
&BEGIN_PROVIDER [ double precision, mo_left_grad_at_fitcusp_radius , (mo_num_8,nucl_num) ]
&BEGIN_PROVIDER [ double precision, mo_left_lapl_at_fitcusp_radius , (mo_num_8,nucl_num) ]

  BEGIN_DOC
  ! Values of the molecular orbitals without S components on atoms
  END_DOC

  implicit none
  integer :: i, j, k, l

  do k = 1, nucl_num
    do j = 1, mo_num
      mo_left_value_at_fitcusp_radius(j,k) = 0.d0
      mo_left_grad_at_fitcusp_radius (j,k) = 0.d0
      mo_left_lapl_at_fitcusp_radius (j,k) = 0.d0
      !DIR$ VECTOR ALIGNED
      do i = 1, ao_num
        mo_left_value_at_fitcusp_radius(j,k) = mo_left_value_at_fitcusp_radius(j,k) + mo_left_coef(i,j) * ao_value_at_fitcusp_radius(i,k)
        mo_left_grad_at_fitcusp_radius (j,k) = mo_left_grad_at_fitcusp_radius (j,k) + mo_left_coef(i,j) * ao_grad_at_fitcusp_radius (i,k)
        mo_left_lapl_at_fitcusp_radius (j,k) = mo_left_lapl_at_fitcusp_radius (j,k) + mo_left_coef(i,j) * ao_lapl_at_fitcusp_radius (i,k)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, nucl_left_fitcusp_param, (4,mo_num,nucl_num) ]

  BEGIN_DOC
  ! Parameters of the splines
  END_DOC

  implicit none
  character*(80)   :: message
  integer          :: i, k, niter
  double precision :: lap_phi, grad_phi, phi, eta
  double precision :: r, Z

  nucl_left_fitcusp_param = 0.d0
  do k = 1, nucl_num
    Z = nucl_charge(k)
    if (Z < 1.d-2) then
      ! Avoid dummy atoms
      cycle
    endif
    R = nucl_fitcusp_radius(k)
    do i = 1, mo_num

      lap_phi  = mo_left_lapl_at_fitcusp_radius (i,k)
      grad_phi = mo_left_grad_at_fitcusp_radius (i,k)
      phi      = mo_left_value_at_fitcusp_radius(i,k)
      eta      = mo_left_value_at_nucl          (i,k)

      nucl_left_fitcusp_param(1,i,k) = -(R*(2.d0*eta*Z-6.d0*grad_phi)+lap_phi*R*R+6.d0*phi)/(2.d0*R*Z-6.d0)
      nucl_left_fitcusp_param(2,i,k) = (lap_phi*R*R*Z-6.d0*grad_phi*R*Z+6.d0*phi*Z+6.d0*eta*Z)/(2.d0*R*Z-6.d0)
      nucl_left_fitcusp_param(3,i,k) = -(R*(-5.d0*grad_phi*Z-1.5d0*lap_phi)+lap_phi*R*R*Z+3.d0*phi*Z+&
          3.d0*eta*Z+6.d0*grad_phi)/(R*R*Z-3.d0*R)
      nucl_left_fitcusp_param(4,i,k) = (R*(-2.d0*grad_phi*Z-lap_phi)+0.5d0*lap_phi*R*R*Z+phi*Z+&
          eta*Z+3.d0*grad_phi)/(R*R*R*Z-3.d0*R*R)

    enddo
  enddo

END_PROVIDER

! ---


