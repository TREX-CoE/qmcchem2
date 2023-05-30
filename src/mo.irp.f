 BEGIN_PROVIDER [ integer, mo_num ]
&BEGIN_PROVIDER [ integer, mo_num_8 ]
  implicit none
  BEGIN_DOC
! Number of Molecular orbitals
  END_DOC
  integer, external              :: mod_align

  mo_num = maxval(present_mos)
  call iinfo(irp_here,'mo_num',mo_num)

  mo_num_8 = mod_align(mo_num)

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, mo_coef_input, (ao_num_8,mo_tot_num) ]

  BEGIN_DOC
  ! Molecular orbital coefficients read from the input file
  END_DOC

  implicit none
  integer           :: i, j
  real, allocatable :: buffer(:,:)

  allocate(buffer(ao_num,mo_tot_num))
  buffer = 0.

  call get_mo_basis_mo_coef(buffer)
  !if(sgn_jast .eq. (-1.d0)) then
  !  call get_bi_ortho_mos_mo_l_coef(buffer)
  !else
  !  call get_mo_basis_mo_coef(buffer)
  !endif

  do i = 1, mo_tot_num
    do j = 1, ao_num
      mo_coef_input(j,i) = buffer(j,i)
    enddo
    call set_order(mo_coef_input(1,i),ao_nucl_sort_idx,ao_num)
    do j = ao_num+1, ao_num_8
      mo_coef_input(j,i) = 0.
    enddo
  enddo

  deallocate(buffer)

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, mo_coef, (ao_num_8,mo_num_8) ]
  implicit none
  BEGIN_DOC
! Molecular orbital coefficients
  END_DOC
  integer                        :: i, j

  do j=1,mo_num
    do i=1,ao_num_8
      mo_coef(i,j) = mo_coef_input(i,j)
    enddo
    !call set_order(mo_coef(1,j),ao_nucl_sort_idx,ao_num)
  enddo
  do j =mo_num+1,mo_num_8
    !DIR$ VECTOR ALIGNED
    do i=1,ao_num_8
      mo_coef(i,j) = 0.
    enddo
  enddo

  ! Input MOs are not needed any more
  FREE mo_coef_input

END_PROVIDER

! ---


BEGIN_PROVIDER [ real, mo_coef_transp, (mo_num_8,ao_num_8) ]
  implicit none
  BEGIN_DOC
! Transpose of the Molecular orbital coefficients
  END_DOC
  call transpose(mo_coef,ao_num_8,mo_coef_transp,mo_num_8,ao_num_8,mo_num_8)

END_PROVIDER


 BEGIN_PROVIDER [ integer, mo_coef_transp_non_zero_idx, (0:mo_num,ao_num) ]
&BEGIN_PROVIDER [ real, mo_coef_transp_sparsity ]
  implicit none
  BEGIN_DOC
! Indices of the non-zero elements of the transpose of the Molecular
! orbital coefficients
  END_DOC
  integer                        :: i, j

  integer                        :: idx
  mo_coef_transp_sparsity = 0.
  do j=1,ao_num
    idx = 0
    do i=1,mo_num
      if (mo_coef_transp(i,j) /= 0.) then
        idx += 1
        mo_coef_transp_non_zero_idx(idx,j) = i
      endif
    enddo
    mo_coef_transp_non_zero_idx(0,j) = idx
    mo_coef_transp_sparsity += float(idx)
  enddo
  mo_coef_transp_sparsity *= 1./(mo_num*ao_num)

END_PROVIDER


BEGIN_PROVIDER [ real, mo_coef_transp_present, (num_present_mos_8,ao_num_8) ]
  implicit none
  BEGIN_DOC
! mo_coef_transp without MOs absent in all determinants
  END_DOC
  integer                        :: i,j,n
  mo_coef_transp_present = 0.
  do i=1,ao_num
    do j=1,num_present_mos
      mo_coef_transp_present(j,i) = mo_coef_transp(present_mos(j),i)
    enddo
  enddo

END_PROVIDER


 BEGIN_PROVIDER [ real, mo_value_transp , (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_grad_transp_x, (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_grad_transp_y, (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_grad_transp_z, (mo_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, mo_lapl_transp  , (mo_num_8,elec_num) ]
  implicit none

  BEGIN_DOC
! Values, gradients, laplacians of the molecular orbitals
!
! Arrays are padded for efficiency
  END_DOC

  integer                        :: i, j, k, l, m

  PROVIDE primitives_reduced

  if (use_qmckl) then
    PROVIDE qmckl_ctx

    do i=1,elec_num
      mo_value_transp (1:num_present_mos,i) = qmckl_mo_vgl(1:num_present_mos,1,i)
      mo_grad_transp_x(1:num_present_mos,i) = qmckl_mo_vgl(1:num_present_mos,2,i)
      mo_grad_transp_y(1:num_present_mos,i) = qmckl_mo_vgl(1:num_present_mos,3,i)
      mo_grad_transp_z(1:num_present_mos,i) = qmckl_mo_vgl(1:num_present_mos,4,i)
      mo_lapl_transp  (1:num_present_mos,i) = qmckl_mo_vgl(1:num_present_mos,5,i)
    end do

  else

    if (do_nucl_fitcusp) then
      PROVIDE nucl_fitcusp_param
      PROVIDE nucl_elec_dist_vec
      PROVIDE nucl_elec_dist_inv
      mo_value_transp  = 0.
      mo_grad_transp_x = 0.
      mo_grad_transp_y = 0.
      mo_grad_transp_z = 0.
      mo_lapl_transp   = 0.
    endif

    do i=1,elec_num
      if (i>1) then
        ao_elec = i
        TOUCH ao_elec
      endif
      if (num_present_mos == mo_num) then
        call sparse_full_mv(mo_coef_transp,mo_num_8,                   &
            ao_value_block(1),ao_num_8,                                &
            ao_grad_block_x(1),                                        &
            ao_grad_block_y(1),                                        &
            ao_grad_block_z(1),                                        &
            ao_lapl_block(1),                                          &
            ao_value_non_zero_idx(0),                                  &
            mo_value_transp(1,i),mo_num_8,                             &
            mo_grad_transp_x(1,i),                                     &
            mo_grad_transp_y(1,i),                                     &
            mo_grad_transp_z(1,i),                                     &
            mo_lapl_transp(1,i),                                       &
            ao_num)

      else
        call sparse_full_mv(mo_coef_transp_present,num_present_mos_8,  &
            ao_value_block(1),ao_num_8,                                &
            ao_grad_block_x(1),                                        &
            ao_grad_block_y(1),                                        &
            ao_grad_block_z(1),                                        &
            ao_lapl_block(1),                                          &
            ao_value_non_zero_idx(0),                                  &
            mo_value_transp(1,i),mo_num_8,                             &
            mo_grad_transp_x(1,i),                                     &
            mo_grad_transp_y(1,i),                                     &
            mo_grad_transp_z(1,i),                                     &
            mo_lapl_transp(1,i),                                       &
            ao_num)

      endif
    enddo  ! i

    ao_elec = 1
    SOFT_TOUCH ao_elec

    if (do_nucl_fitcusp) then
      real                           :: r, r_inv, c
      do i=1,elec_num
        do k=1,nucl_num
          r = nucl_elec_dist(k,i)
          if (r > nucl_fitcusp_radius(k)) then
            cycle
          endif
          r_inv = nucl_elec_dist_inv(k,i)
          do j=1,num_present_mos
            mo_value_transp(j,i) = mo_value_transp(j,i) +              &
                     nucl_fitcusp_param(1,j,k) +                       &
                r * (nucl_fitcusp_param(2,j,k) +                       &
                r * (nucl_fitcusp_param(3,j,k) +                       &
                r *  nucl_fitcusp_param(4,j,k) ))
            mo_lapl_transp(j,i) = mo_lapl_transp(j,i) +                &
                nucl_fitcusp_param(2,j,k)*(r_inv+r_inv) +              &
                      6.*nucl_fitcusp_param(3,j,k) +                   &
                 r * 12.*nucl_fitcusp_param(4,j,k)
            c = r_inv * (nucl_fitcusp_param(2,j,k) +                   &
                 r * (2.*nucl_fitcusp_param(3,j,k) +                   &
                 r *  3.*nucl_fitcusp_param(4,j,k) ))
            mo_grad_transp_x(j,i) = mo_grad_transp_x(j,i) + nucl_elec_dist_vec(1,k,i)*c
            mo_grad_transp_y(j,i) = mo_grad_transp_y(j,i) + nucl_elec_dist_vec(2,k,i)*c
            mo_grad_transp_z(j,i) = mo_grad_transp_z(j,i) + nucl_elec_dist_vec(3,k,i)*c
          enddo
          ! It is safe to exit here because a core electron is close to only one nucleus
          exit
        enddo ! k
      enddo ! i

    endif

  endif

  if (do_prepare) then
    real                           :: lambda, t
    ! Scale off-diagonal elements
    t = prepare_walkers_t
    do j=1,elec_alpha_num
      do i=1,num_present_mos
        if (i /= j) then
          mo_value_transp(i,j) *= t
          mo_grad_transp_x(i,j) *= t
          mo_grad_transp_y(i,j) *= t
          mo_grad_transp_z(i,j) *= t
          mo_lapl_transp(i,j) *= t
        endif
      enddo
    enddo
    do j=1,elec_beta_num
      do i=1,num_present_mos
        if (i /= j) then
          mo_value_transp(i,j+elec_alpha_num) *= t
          mo_grad_transp_x(i,j+elec_alpha_num) *= t
          mo_grad_transp_y(i,j+elec_alpha_num) *= t
          mo_grad_transp_z(i,j+elec_alpha_num) *= t
          mo_lapl_transp(i,j+elec_alpha_num) *= t
        endif
      enddo
    enddo
  endif

  if (num_present_mos < mo_num) then
    do i=1,elec_num
        do j=num_present_mos,1,-1
          mo_value_transp (present_mos(j),i) = mo_value_transp (j,i)
          mo_grad_transp_x(present_mos(j),i) = mo_grad_transp_x(j,i)
          mo_grad_transp_y(present_mos(j),i) = mo_grad_transp_y(j,i)
          mo_grad_transp_z(present_mos(j),i) = mo_grad_transp_z(j,i)
          mo_lapl_transp  (present_mos(j),i) = mo_lapl_transp  (j,i)
          if (present_mos(j) == j) then
            exit
          endif
        enddo
    enddo
  endif


END_PROVIDER

BEGIN_PROVIDER [ real, mo_value, (elec_num_8,mo_num) ]
  implicit none
  BEGIN_DOC
! Values of the molecular orbitals
  END_DOC

  integer                        :: i,j
  integer, save                  :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    PROVIDE primitives_reduced
    !DIR$ VECTOR ALIGNED
    mo_value = 0.
  endif
  call transpose(mo_value_transp(1,1),mo_num_8,mo_value,elec_num_8,mo_num,elec_num)

END_PROVIDER


 BEGIN_PROVIDER [ double precision, mo_grad_x, (elec_num_8,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_grad_y, (elec_num_8,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_grad_z, (elec_num_8,mo_num) ]
  implicit none

  BEGIN_DOC
! Gradients of the molecular orbitals
  END_DOC

  integer                        :: i,j
  integer, save                  :: ifirst = 0

  if (ifirst == 0) then
    !DIR$ VECTOR ALIGNED
    mo_grad_x = 0.d0
    !DIR$ VECTOR ALIGNED
    mo_grad_y = 0.d0
    !DIR$ VECTOR ALIGNED
    mo_grad_z = 0.d0
    ifirst = 1
    PROVIDE primitives_reduced
  endif
  ! Transpose x last for cache efficiency
  call transpose_to_dp(mo_grad_transp_y(1,1),mo_num_8,mo_grad_y(1,1),elec_num_8,mo_num,elec_num)
  call transpose_to_dp(mo_grad_transp_z(1,1),mo_num_8,mo_grad_z(1,1),elec_num_8,mo_num,elec_num)
  call transpose_to_dp(mo_grad_transp_x(1,1),mo_num_8,mo_grad_x(1,1),elec_num_8,mo_num,elec_num)


END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_lapl, (elec_num_8,mo_num) ]
  implicit none
  BEGIN_DOC
! Laplacians of the molecular orbitals
  END_DOC

  integer                        :: i,j
  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    PROVIDE primitives_reduced
    !DIR$ VECTOR ALIGNED
    mo_lapl = 0.d0
  endif
  call transpose_to_dp(mo_lapl_transp(1,1),mo_num_8,mo_lapl,elec_num_8,mo_num,elec_num)

END_PROVIDER


 BEGIN_PROVIDER [ double precision, mo_grad_lapl_alpha, (4,elec_alpha_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_grad_lapl_beta , (4,elec_alpha_num+1:elec_num,mo_num) ]
 implicit none
 BEGIN_DOC
! Gradients and laplacian
 END_DOC
 integer :: i,j
 do j=1,mo_num
    do i=1,elec_alpha_num
      mo_grad_lapl_alpha(1,i,j) = mo_grad_transp_x(j,i)
      mo_grad_lapl_alpha(2,i,j) = mo_grad_transp_y(j,i)
      mo_grad_lapl_alpha(3,i,j) = mo_grad_transp_z(j,i)
      mo_grad_lapl_alpha(4,i,j) = mo_lapl_transp  (j,i)
    enddo
  enddo
  do j=1,mo_num
    do i=elec_alpha_num+1,elec_num
      mo_grad_lapl_beta(1,i,j) = mo_grad_transp_x(j,i)
      mo_grad_lapl_beta(2,i,j) = mo_grad_transp_y(j,i)
      mo_grad_lapl_beta(3,i,j) = mo_grad_transp_z(j,i)
      mo_grad_lapl_beta(4,i,j) = mo_lapl_transp  (j,i)
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_grad_lapl_transp, (4,mo_num,elec_num) ]
 implicit none
 BEGIN_DOC
! Gradients and laplacian
 END_DOC
 integer :: i,j
 do i=1,elec_num
   do j=1,mo_num
     mo_grad_lapl_transp(1,j,i) = mo_grad_transp_x(j,i)
     mo_grad_lapl_transp(2,j,i) = mo_grad_transp_y(j,i)
     mo_grad_lapl_transp(3,j,i) = mo_grad_transp_z(j,i)
     mo_grad_lapl_transp(4,j,i) = mo_lapl_transp  (j,i)
   enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ real, prepare_walkers_t ]
  implicit none
  BEGIN_DOC
! prepare_walkers_t : scaling of the off-diagonal elements
! of the mo_value matrix
  END_DOC
  prepare_walkers_t = 1.
END_PROVIDER


BEGIN_PROVIDER [ integer, mo_tot_num ]

  BEGIN_DOC
! Total number of MOs in the EZFIO file
  END_DOC

  mo_tot_num = -1
  call get_mo_basis_mo_num(mo_tot_num)
  if (mo_tot_num <= 0) then
    call abrt(irp_here,'Total number of MOs can''t be <0')
  endif
  call iinfo(irp_here,'mo_tot_num',mo_tot_num)

END_PROVIDER

! ---

BEGIN_PROVIDER [ integer, n_oo ]

  !if(mod(elec_num, 2) .eq. 0) then
  if(elec_alpha_num .eq. elec_beta_num) then
    n_oo = elec_alpha_num
  else
    n_oo = max(elec_alpha_num, elec_beta_num)
  endif

END_PROVIDER


!-----------------
! Fit cusp
!-----------------


BEGIN_PROVIDER [ double precision , mo_value_at_nucl, (mo_num_8,nucl_num) ]
  implicit none
  BEGIN_DOC
! Values of the molecular orbitals at the nucleus without the
! S components of the current nucleus
  END_DOC
  integer                        :: i, j, k, l
  real                           :: ao_value_at_nucl_no_S(ao_num)

  PROVIDE mo_fitcusp_normalization_before
  do k=1,nucl_num
    point(1) = nucl_coord(k,1)
    point(2) = nucl_coord(k,2)
    point(3) = nucl_coord(k,3)
    TOUCH point

    PROVIDE ao_value_p

    do i=1,ao_num
      if (ao_nucl(i) /= k) then
        ao_value_at_nucl_no_S(i) = ao_value_p(i)
      else
        ao_value_at_nucl_no_S(i) = 0.
      endif
    enddo

    integer :: jj
    do jj=1,num_present_mos
      j = present_mos(jj)
      mo_value_at_nucl(j,k) = 0.
      !DIR$ VECTOR ALIGNED
      do i=1,ao_num
        mo_value_at_nucl(j,k) = mo_value_at_nucl(j,k) + mo_coef(i,j)*ao_value_at_nucl_no_S(i)
      enddo
    enddo

  enddo
  FREE ao_value_p ao_grad_p ao_lapl_p ao_axis_grad_p ao_oned_grad_p ao_oned_prim_grad_p
  FREE ao_oned_lapl_p ao_axis_lapl_p ao_oned_prim_lapl_p ao_oned_p ao_oned_prim_p
  FREE ao_axis_p ao_axis_power_p
  SOFT_TOUCH point

END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_value_at_fitcusp_radius, (ao_num_8,nucl_num) ]
&BEGIN_PROVIDER [ double precision, ao_grad_at_fitcusp_radius, (ao_num_8,nucl_num) ]
&BEGIN_PROVIDER [ double precision, ao_lapl_at_fitcusp_radius, (ao_num_8,nucl_num) ]
  implicit none
  BEGIN_DOC
! Values of the atomic orbitals with only S components on atoms
  END_DOC

  integer                        :: i, j, k

  do k=1,nucl_num
    point(1) = nucl_coord(k,1)
    point(2) = nucl_coord(k,2)
    point(3) = nucl_coord(k,3)+ nucl_fitcusp_radius(k)
    TOUCH point

    do j=1,ao_num
      ao_value_at_fitcusp_radius(j,k) = ao_value_p(j)
      ao_grad_at_fitcusp_radius(j,k) = ao_grad_p(j,3)
      ao_lapl_at_fitcusp_radius(j,k) = ao_lapl_p(j)
      if ( (ao_nucl(j) /= k).or.(ao_power(j,4) >0) ) then
        ao_value_at_fitcusp_radius(j,k) = 0.
        ao_grad_at_fitcusp_radius(j,k) = 0.
        ao_lapl_at_fitcusp_radius(j,k) = 0.
      endif
    enddo
  enddo
  FREE ao_value_p ao_grad_p ao_lapl_p ao_axis_grad_p ao_oned_grad_p ao_oned_prim_grad_p
  FREE ao_oned_lapl_p ao_axis_lapl_p ao_oned_prim_lapl_p ao_oned_p ao_oned_prim_p
  FREE ao_axis_p ao_axis_power_p
  SOFT_TOUCH point

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_fitcusp_normalization_before, (mo_tot_num) ]
 implicit none
 BEGIN_DOC
 ! Renormalization factor of MOs due to cusp fitting
 END_DOC
 include 'constants.F'
 integer :: i,j,k,l
 double precision :: dr, r, f, t
 integer, save :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    mo_fitcusp_normalization_before = 0.d0
    do k=1,nucl_num
      dr = nucl_fitcusp_radius(k)*1.d-2
      point(1) = nucl_coord(k,1)
      point(2) = nucl_coord(k,2)
      point(3) = nucl_coord(k,3)-dr
      do l=1,101
        r = point(3) + dr
        point(3) = r
        TOUCH point
        f = dfour_pi*r*r*dr
        do i=1,mo_tot_num
          mo_fitcusp_normalization_before(i) += f*mo_value_p(i)**2
        enddo
      enddo
    enddo
  endif

END_PROVIDER


BEGIN_PROVIDER [ double precision, mo_fitcusp_normalization_after, (mo_tot_num) ]
 implicit none
 BEGIN_DOC
 ! Renormalization factor of MOs due to cusp fitting
 END_DOC
 include 'constants.F'
 integer :: i,j,k,l
 double precision :: dr, r, f, t, t2
 integer, save :: ifirst = 0

  PROVIDE primitives_reduced
  if (ifirst == 0) then
    ifirst = 1
    mo_fitcusp_normalization_after = 0.d0
    do k=1,nucl_num
      dr = nucl_fitcusp_radius(k)*1.d-2
      point(1) = nucl_coord(k,1)
      point(2) = nucl_coord(k,2)
      point(3) = nucl_coord(k,3)- dr
      do l=1,101
        point(3) = point(3)+ dr
        TOUCH point nucl_fitcusp_param primitives_reduced mo_coef
        r = point(3)
        f = dfour_pi*r*r*dr
        do i=1,mo_num
          t = 0.d0
          do j=1,ao_num
            if ( (ao_nucl(j) /= k).or.(ao_power(j,4) > 0) ) then
              t = t + mo_coef(j,i) * ao_value_p(j)
            endif
          enddo
          t = t + nucl_fitcusp_param(1,i,k) +                        &
              r * (nucl_fitcusp_param(2,i,k) +                       &
              r * (nucl_fitcusp_param(3,i,k) +                       &
              r *  nucl_fitcusp_param(4,i,k) ))
          mo_fitcusp_normalization_after(i) += t*t*f
        enddo
      enddo
    enddo
  endif

END_PROVIDER


 BEGIN_PROVIDER [ double precision, mo_value_at_fitcusp_radius, (num_present_mos_8,nucl_num) ]
&BEGIN_PROVIDER [ double precision, mo_grad_at_fitcusp_radius, (num_present_mos_8,nucl_num) ]
&BEGIN_PROVIDER [ double precision, mo_lapl_at_fitcusp_radius, (num_present_mos_8,nucl_num) ]
  implicit none
  BEGIN_DOC
! Values of the molecular orbitals without S components on atoms.
  END_DOC
  integer                        :: i, j, k, l

   do k=1,nucl_num
     do j=1,num_present_mos
       mo_value_at_fitcusp_radius(j,k) = 0.d0
       mo_grad_at_fitcusp_radius(j,k) = 0.d0
       mo_lapl_at_fitcusp_radius(j,k) = 0.d0
     enddo
     do i=1,ao_num
       !DIR$ VECTOR ALIGNED
       do j=1,num_present_mos
         mo_value_at_fitcusp_radius(j,k) = mo_value_at_fitcusp_radius(j,k) + mo_coef_transp_present(j,i)*ao_value_at_fitcusp_radius(i,k)
         mo_grad_at_fitcusp_radius (j,k) = mo_grad_at_fitcusp_radius (j,k) + mo_coef_transp_present(j,i)*ao_grad_at_fitcusp_radius (i,k)
         mo_lapl_at_fitcusp_radius (j,k) = mo_lapl_at_fitcusp_radius (j,k) + mo_coef_transp_present(j,i)*ao_lapl_at_fitcusp_radius (i,k)
       enddo
     enddo
   enddo

END_PROVIDER


BEGIN_PROVIDER [ real, nucl_fitcusp_param, (4,mo_num,nucl_num) ]
  implicit none
  BEGIN_DOC
! Parameters of the splines
  END_DOC
  integer                        :: i,k, niter
  character*(80)                 :: message

  nucl_fitcusp_param = 0.d0
  do k=1,nucl_num
    double precision               :: r, Z
    Z = nucl_charge(k)
    if (Z < 1.d-2) then
      ! Avoid dummy atoms
      cycle
    endif
    R = nucl_fitcusp_radius(k)
    do i=1,mo_num
      double precision               :: lap_phi, grad_phi, phi, eta
      lap_phi = mo_lapl_at_fitcusp_radius(i,k)
      grad_phi = mo_grad_at_fitcusp_radius(i,k)
      phi = mo_value_at_fitcusp_radius(i,k)
      eta = mo_value_at_nucl(i,k)

      nucl_fitcusp_param(1,i,k) = -(R*(2.d0*eta*Z-6.d0*grad_phi)+lap_phi*R*R+6.d0*phi)/(2.d0*R*Z-6.d0)
      nucl_fitcusp_param(2,i,k) = (lap_phi*R*R*Z-6.d0*grad_phi*R*Z+6.d0*phi*Z+6.d0*eta*Z)/(2.d0*R*Z-6.d0)
      nucl_fitcusp_param(3,i,k) = -(R*(-5.d0*grad_phi*Z-1.5d0*lap_phi)+lap_phi*R*R*Z+3.d0*phi*Z+&
          3.d0*eta*Z+6.d0*grad_phi)/(R*R*Z-3.d0*R)
      nucl_fitcusp_param(4,i,k) = (R*(-2.d0*grad_phi*Z-lap_phi)+0.5d0*lap_phi*R*R*Z+phi*Z+&
          eta*Z+3.d0*grad_phi)/(R*R*R*Z-3.d0*R*R)

    enddo
  enddo

END_PROVIDER



subroutine sparse_full_mv(A,LDA,                                     &
      B1,LDB,                                                        &
      B2, B3, B4, B5, indices,                                       &
      C1,LDC,C2,C3,C4,C5,an)
  implicit none
  BEGIN_DOC
! Performs a vectorized product between a dense matrix (the MO coefficients
! matrix) and 5 sparse vectors (the value, gradients and laplacian of the AOs).
  END_DOC
  integer, intent(in)            :: an,LDA,LDB,LDC
  integer, intent(in)            :: indices(0:LDB)
  real, intent(in)               :: A(LDA,an)
  real, intent(in)               :: B1(LDB)
  real, intent(in)               :: B2(LDB)
  real, intent(in)               :: B3(LDB)
  real, intent(in)               :: B4(LDB)
  real, intent(in)               :: B5(LDB)
  real, intent(out)              :: C1(LDC)
  real, intent(out)              :: C2(LDC)
  real, intent(out)              :: C3(LDC)
  real, intent(out)              :: C4(LDC)
  real, intent(out)              :: C5(LDC)
  !DIR$ ASSUME_ALIGNED A  : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED B1 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED B2 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED B3 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED B4 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED B5 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED C1 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED C2 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED C3 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED C4 : $IRP_ALIGN
  !DIR$ ASSUME_ALIGNED C5 : $IRP_ALIGN

  integer                        :: kao, kmax, kmax2, kmax3
  integer                        :: i,j,k
  integer                        :: k_vec(8)
  !DIR$ ATTRIBUTES ALIGN: $IRP_ALIGN :: k_vec
  real                           :: d11, d12, d13, d14, d15
  real                           :: d21, d22, d23, d24, d25
  real                           :: d31, d32, d33, d34, d35
  real                           :: d41, d42, d43, d44, d45

  ! LDC and LDA have to be factors of simd_sp

!  IRP_IF NO_PREFETCH
!  IRP_ELSE
!    call MM_PREFETCH (A(1,indices(1)),3)
!    call MM_PREFETCH (A(1,indices(2)),3)
!    call MM_PREFETCH (A(1,indices(3)),3)
!    call MM_PREFETCH (A(1,indices(4)),3)
!  IRP_ENDIF

  !OMP$ SIMD
  do j=1,LDC
    C1(j) = 0.
    C2(j) = 0.
    C3(j) = 0.
    C4(j) = 0.
    C5(j) = 0.
  enddo

  kmax2 = shiftr(indices(0),2)
  kmax2 = shiftl(kmax2,2)
  kmax3 = indices(0)

  do kao=1,kmax2,4
    k_vec(1) = indices(kao  )
    k_vec(2) = indices(kao+1)
    k_vec(3) = indices(kao+2)
    k_vec(4) = indices(kao+3)

    d11 = B1(kao  )
    d21 = B1(kao+1)
    d31 = B1(kao+2)
    d41 = B1(kao+3)

    d12 = B2(kao  )
    d22 = B2(kao+1)
    d32 = B2(kao+2)
    d42 = B2(kao+3)

    d13 = B3(kao  )
    d23 = B3(kao+1)
    d33 = B3(kao+2)
    d43 = B3(kao+3)

    d14 = B4(kao  )
    d24 = B4(kao+1)
    d34 = B4(kao+2)
    d44 = B4(kao+3)

    d15 = B5(kao  )
    d25 = B5(kao+1)
    d35 = B5(kao+2)
    d45 = B5(kao+3)

    do k=0,LDA-1,$IRP_ALIGN/4
      !DIR$ VECTOR ALIGNED
      !OMP$ SIMD FIRSTPRIVATE(d11,d21,d31,d41)
      do j=1,$IRP_ALIGN/4
!        IRP_IF NO_PREFETCH
!        IRP_ELSE
!          call MM_PREFETCH (A(j+k,indices(kao+4)),3)
!          call MM_PREFETCH (A(j+k,indices(kao+5)),3)
!          call MM_PREFETCH (A(j+k,indices(kao+6)),3)
!          call MM_PREFETCH (A(j+k,indices(kao+7)),3)
!        IRP_ENDIF
        C1(j+k) = C1(j+k) + A(j+k,k_vec(1))*d11 + A(j+k,k_vec(2))*d21&
            +   A(j+k,k_vec(3))*d31 + A(j+k,k_vec(4))*d41
      enddo

      !DIR$ VECTOR ALIGNED
      !OMP$ SIMD FIRSTPRIVATE(d12,d22,d32,d42,d13,d23,d33,d43)
      do j=1,$IRP_ALIGN/4
        C2(j+k) = C2(j+k) + A(j+k,k_vec(1))*d12 + A(j+k,k_vec(2))*d22&
            +   A(j+k,k_vec(3))*d32 + A(j+k,k_vec(4))*d42
        C3(j+k) = C3(j+k) + A(j+k,k_vec(1))*d13 + A(j+k,k_vec(2))*d23&
            +   A(j+k,k_vec(3))*d33 + A(j+k,k_vec(4))*d43
      enddo

      !DIR$ VECTOR ALIGNED
      !OMP$ SIMD FIRSTPRIVATE(d14,d24,d34,d44,d15,d25,d35,d45)
      do j=1,$IRP_ALIGN/4
        C4(j+k) = C4(j+k) + A(j+k,k_vec(1))*d14 + A(j+k,k_vec(2))*d24&
            +   A(j+k,k_vec(3))*d34 + A(j+k,k_vec(4))*d44
        C5(j+k) = C5(j+k) + A(j+k,k_vec(1))*d15 + A(j+k,k_vec(2))*d25&
            +   A(j+k,k_vec(3))*d35 + A(j+k,k_vec(4))*d45
      enddo
    enddo

  enddo

  do kao = kmax2+1, kmax3
    k_vec(1) = indices(kao)
    d11 = B1(kao)
    d12 = B2(kao)
    d13 = B3(kao)
    d14 = B4(kao)
    d15 = B5(kao)
    !DIR$ VECTOR ALIGNED
    do k=0,LDA-1,$IRP_ALIGN/4
      !DIR$ VECTOR ALIGNED
      !OMP$ SIMD FIRSTPRIVATE(d11,d12,d13,d14,d15)
      do j=1,$IRP_ALIGN/4
        C1(j+k) = C1(j+k) + A(j+k,k_vec(1))*d11
        C2(j+k) = C2(j+k) + A(j+k,k_vec(1))*d12
        C3(j+k) = C3(j+k) + A(j+k,k_vec(1))*d13
        C4(j+k) = C4(j+k) + A(j+k,k_vec(1))*d14
        C5(j+k) = C5(j+k) + A(j+k,k_vec(1))*d15
      enddo
    enddo
  enddo

end



