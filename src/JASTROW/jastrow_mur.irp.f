! Mu(r) Jastrow
! --------------

! ---

 BEGIN_PROVIDER [double precision, jast_Mur_value    ]
&BEGIN_PROVIDER [double precision, jast_Mur_value_inv]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Mur_value(i)
  enddo

  jast_Mur_value     = dexp(argexpo)
  jast_Mur_value_inv = 1.d0 / jast_Mur_value

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Mur_value , (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mur_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mur_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mur_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mur_lapl  , (elec_num_8)]

  implicit none
  integer          :: i, j
  double precision :: r1(3), r2(3), je_val, je_der(3), je_lap
  double precision :: dx, dy, dz, r12

  do i = 1, elec_num

    r1(1) = elec_coord(i,1)
    r1(2) = elec_coord(i,2)
    r1(3) = elec_coord(i,3)

    jast_elec_Mur_value (i) = 0.d0
    jast_elec_Mur_grad_x(i) = 0.d0
    jast_elec_Mur_grad_y(i) = 0.d0
    jast_elec_Mur_grad_z(i) = 0.d0
    jast_elec_Mur_lapl  (i) = 0.d0

    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num

      if(i==j) cycle

      r2(1) = elec_coord(j,1)
      r2(2) = elec_coord(j,2)
      r2(3) = elec_coord(j,3)

      !dx  = r1(1) - r2(1)
      !dy  = r1(2) - r2(2)
      !dz  = r1(3) - r2(3)
      !r12 = dsqrt(dx*dx + dy*dy + dz*dz)
      !je_val    = r12
      !je_der(1) = dx / r12
      !je_der(2) = dy / r12
      !je_der(3) = dz / r12
      !je_lap    = 0.d0
      call j_elec_mur(r1, r2, je_val, je_der, je_lap)

      jast_elec_Mur_value (i) = jast_elec_Mur_value (i) + 0.5d0 * je_val
      jast_elec_Mur_grad_x(i) = jast_elec_Mur_grad_x(i) + je_der(1)
      jast_elec_Mur_grad_y(i) = jast_elec_Mur_grad_y(i) + je_der(2)
      jast_elec_Mur_grad_z(i) = jast_elec_Mur_grad_z(i) + je_der(3)
      jast_elec_Mur_lapl  (i) = jast_elec_Mur_lapl  (i) + je_lap
    enddo
  enddo

END_PROVIDER

! ---

subroutine j_elec_mur(r1, r2, je_val, je_der, je_lap)

  BEGIN_DOC  
  !
  ! J(r1,r2) = 0.5 [ r12 (1-erf(mu(r1,r2) r12)) - exp(-(mu(r1,r2) r12)**2) / (pi**0.5 mu(r1,r2)) ] 
  !
  END_DOC

  include '../constants.F'

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: je_val, je_der(3), je_lap
  double precision              :: mur_val, mur_der(3), mur_lap
  double precision              :: r12, dx, dy, dz
  double precision              :: tmp, mur_tmp, mur_der_square

  dx  = r1(1) - r2(1)
  dy  = r1(2) - r2(2) 
  dz  = r1(3) - r2(3)
  r12 = dsqrt(dx*dx + dy*dy + dz*dz) + 1d-12

  call mur_val_der_lap_r1(r1, r2, mur_val, mur_der, mur_lap)

  mur_tmp = mur_val * r12
  je_val  = 0.5d0 * (r12 * (1.d0 - derf(mur_tmp)) - dexp(-mur_tmp*mur_tmp)/(dsqpi*mur_val))

  tmp       = 0.5d0 * (1.d0 - derf(mur_tmp)) / r12 
  je_der(1) = tmp * dx
  je_der(2) = tmp * dy
  je_der(3) = tmp * dz
  je_lap    = 2.d0 * tmp - (mur_val + 2.d0*(mur_der(1)*dx + mur_der(2)*dy + mur_der(3)*dz)) * dexp(-mur_tmp*mur_tmp) / dsqpi

  ! NB
  ! no need to dabs because mur_val should be positive
  if(mur_val .gt. 1d-7) then
    tmp            = 0.5d0 * dexp(-mur_tmp*mur_tmp) / (dsqpi*mur_val*mur_val)
    mur_der_square = mur_der(1) * mur_der(1) + mur_der(2) * mur_der(2) + mur_der(3) * mur_der(3)
    je_der(1)      = je_der(1) + tmp * mur_der(1)
    je_der(2)      = je_der(2) + tmp * mur_der(2)
    je_der(3)      = je_der(3) + tmp * mur_der(3)
    je_lap         = je_lap    + tmp * (mur_lap - 2.d0 * mur_val * mur_der_square * (r12*r12 + 1.d0/(mur_val*mur_val)))
  endif

  return
end subroutine j_elec_mur

! ---

subroutine mur_val_der_lap_r1(r1, r2, mur_val, mur_der, mur_lap)

  BEGIN_DOC  
  !
  ! mu(r1,r2)
  ! grad_{r1} mu(r1, r2)
  ! lap_{r1}  mu(r1, r2)
  !
  END_DOC

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: mur_val, mur_der(3), mur_lap
  double precision              :: r(3)
  double precision              :: rho_val, rho_der(3), rho_lap, rho_der_square
  double precision              :: tmp, tmp1, tmp2
  double precision              :: r12, dx, dy, dz

  !dx         = r1(1) - r2(1)
  !dy         = r1(2) - r2(2) 
  !dz         = r1(3) - r2(3)
  !r12        = dsqrt(dx*dx + dy*dy + dz*dz) + 1d-12
  !mur_val    = r12
  !mur_der(1) = dx / r12
  !mur_der(2) = dy / r12
  !mur_der(3) = dz / r12
  !mur_lap    = 2.d0 / r12
  !return

  PROVIDE j1b_type

  if(j1b_type .eq. 200) then

    ! mu[rho(r)] = mu_r_ct \sqrt(rho(r)) + mu_erf \exp(-rho(r))

    PROVIDE mu_r_ct mu_erf

    r(1) = 0.5d0 * (r1(1) + r2(1))
    r(2) = 0.5d0 * (r1(2) + r2(2))
    r(3) = 0.5d0 * (r1(3) + r2(3))

    call rho_hf_val_der_lap(r, rho_val, rho_der, rho_lap)

    mur_val = mu_r_ct * dsqrt(rho_val) + mu_erf * dexp(-rho_val)
  
    tmp        = -0.5d0 * mu_erf * dexp(-rho_val)
    mur_der(1) = tmp * rho_der(1)
    mur_der(2) = tmp * rho_der(2)
    mur_der(3) = tmp * rho_der(3)

    rho_der_square  = rho_der(1) * rho_der(1) + rho_der(2) * rho_der(2) + rho_der(3) * rho_der(3)
    mur_lap         = 0.5d0 * tmp * (rho_lap - rho_der_square)

    ! NB
    ! no need to dabs because rho_val should be positive
    if(rho_val .lt. 1d-8) return
    tmp        = 0.25d0 * mu_r_ct / dsqrt(rho_val)
    mur_der(1) = mur_der(1) + tmp * rho_der(1)
    mur_der(2) = mur_der(2) + tmp * rho_der(2)
    mur_der(3) = mur_der(3) + tmp * rho_der(3)
    mur_lap    = mur_lap + 0.5d0 * tmp * (rho_lap - 0.5d0 * rho_der_square / rho_val)

  elseif(j1b_type .eq. 201) then

    ! mu[rho(r)] = mu_r_ct rho(r) + mu_erf \exp(-rho(r))

    PROVIDE mu_r_ct mu_erf

    r(1) = 0.5d0 * (r1(1) + r2(1))
    r(2) = 0.5d0 * (r1(2) + r2(2))
    r(3) = 0.5d0 * (r1(3) + r2(3))

    call rho_hf_val_der_lap(r, rho_val, rho_der, rho_lap)

    mur_val = mu_r_ct * rho_val + mu_erf * dexp(-rho_val)
  
    tmp1       = mu_erf * dexp(-rho_val)
    tmp2       = 0.5d0 * (mu_r_ct - tmp1)
    mur_der(1) = tmp2 * rho_der(1)
    mur_der(2) = tmp2 * rho_der(2)
    mur_der(3) = tmp2 * rho_der(3)

    rho_der_square  = rho_der(1) * rho_der(1) + rho_der(2) * rho_der(2) + rho_der(3) * rho_der(3)
    mur_lap         = 0.25d0 * (mu_r_ct * rho_lap + tmp1 * (rho_der_square - rho_lap))

  else
 
    print*, 'j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end subroutine mur_val_der_lap_r1

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Mur_grad_x_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mur_grad_y_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mur_grad_z_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mur_lapl_num  , (elec_num_8)]

  implicit none
  integer          :: i, j
  double precision :: r1(3), r2(3), je_der(3), je_lap
  double precision :: eps, tmp_der, tmp_lap, je_p, je_m, je_0, tmp

  eps     = 1d-3
  tmp_der = 0.5d0 /  eps
  tmp_lap = 1.0d0 / (eps * eps)

  do i = 1, elec_num

    r1(1) = elec_coord(i,1)
    r1(2) = elec_coord(i,2)
    r1(3) = elec_coord(i,3)

    jast_elec_Mur_grad_x_num(i) = 0.d0
    jast_elec_Mur_grad_y_num(i) = 0.d0
    jast_elec_Mur_grad_z_num(i) = 0.d0
    jast_elec_Mur_lapl_num  (i) = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      r2(1) = elec_coord(j,1)
      r2(2) = elec_coord(j,2)
      r2(3) = elec_coord(j,3)
      call j_elec_mur(r1, r2, je_0, je_der, je_lap)

      r1(1) += eps
      call j_elec_mur(r1, r2, je_p, je_der, je_lap)
      r1(1) -= 2.d0 * eps
      call j_elec_mur(r1, r2, je_m, je_der, je_lap)
      r1(1) += eps
      jast_elec_Mur_grad_x_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Mur_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      r1(2) += eps
      call j_elec_mur(r1, r2, je_p, je_der, je_lap)
      r1(2) -= 2.d0 * eps
      call j_elec_mur(r1, r2, je_m, je_der, je_lap)
      r1(2) += eps
      jast_elec_Mur_grad_y_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Mur_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      r1(3) += eps
      call j_elec_mur(r1, r2, je_p, je_der, je_lap)
      r1(3) -= 2.d0 * eps
      call j_elec_mur(r1, r2, je_m, je_der, je_lap)
      r1(3) += eps
      jast_elec_Mur_grad_z_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Mur_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

    enddo
  enddo

END_PROVIDER

