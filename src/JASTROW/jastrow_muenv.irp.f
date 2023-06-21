! Mu Jastrow x envelop
! ---------------------

! ---

 BEGIN_PROVIDER [double precision, jast_Muenv_value    ]
&BEGIN_PROVIDER [double precision, jast_Muenv_value_inv]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Muenv_value(i)
  enddo

  jast_Muenv_value     = dexp(argexpo)
  jast_Muenv_value_inv = 1.d0 / jast_Muenv_value

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, jast_elec_Muenv_value, (elec_num_8)]

  BEGIN_DOC  
  !
  ! J(i) = 0.5  \sum_{j!=i} 0.5 [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ]
  !      x v(riA) x v(rjA)
  !
  !      = 0.25 \sum_{j!=i} [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ]
  !      x v(riA) x v(rjA)
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: mu_rij, mu_pi
  double precision :: rij, tmp_ij

  mu_pi = 1.d0 / (dsqpi * mu_erf)

  do i = 1, elec_num

    tmp_ij = 0.d0
    !DIR$ LOOP COUNT(100)
    do j = 1, elec_num

      if(j==i) cycle

      rij    = elec_dist(j,i)
      mu_rij = mu_erf * rij

      tmp_ij += (rij * (1.d0 - derf(mu_rij)) - mu_pi * dexp(-mu_rij*mu_rij)) * vi_1b_env(j)
    enddo

    jast_elec_Muenv_value(i) = 0.25d0 * tmp_ij * vi_1b_env(i)
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Muenv_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Muenv_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Muenv_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Muenv_lapl  , (elec_num_8)]

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: mu_div_sqrtpi, mu_sqrtpi_inv, rij, mu_rij
  double precision :: tmp0_ij, tmp1_ij, tmp2_ij, tmp3_ij
  double precision :: vj_lapl_uij, vj_derivx_uij, vj_derivy_uij, vj_derivz_uij, vj_uij

  mu_div_sqrtpi = mu_erf / dsqpi
  mu_sqrtpi_inv = 1.d0 / ( dsqpi * mu_erf )

  do i = 1, elec_num

    vj_uij        = 0.d0
    vj_derivx_uij = 0.d0 
    vj_derivy_uij = 0.d0 
    vj_derivz_uij = 0.d0 
    vj_lapl_uij   = 0.d0

    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num

      if(i==j) cycle

      rij     = elec_dist(j,i)
      mu_rij  = mu_erf * rij
      tmp0_ij = dexp(-mu_rij * mu_rij)
      tmp1_ij = 1.d0 - derf(mu_rij)
      tmp2_ij = tmp1_ij * elec_dist_inv(j,i)
      tmp3_ij = -0.5d0 * tmp2_ij * vi_1b_env(j)

      vj_uij        += 0.5d0 * (rij * tmp1_ij - mu_sqrtpi_inv * tmp0_ij) * vi_1b_env(j)
      vj_derivx_uij += tmp3_ij * elec_dist_vec_x(j,i)
      vj_derivy_uij += tmp3_ij * elec_dist_vec_y(j,i)
      vj_derivz_uij += tmp3_ij * elec_dist_vec_z(j,i)
      vj_lapl_uij   += (tmp2_ij - mu_div_sqrtpi * tmp0_ij) * vi_1b_env(j)
    enddo

    jast_elec_Muenv_grad_x(i) = vj_derivx_uij * vi_1b_env(i) + vj_uij * deriv_vi_x_env(i)
    jast_elec_Muenv_grad_y(i) = vj_derivy_uij * vi_1b_env(i) + vj_uij * deriv_vi_y_env(i)
    jast_elec_Muenv_grad_z(i) = vj_derivz_uij * vi_1b_env(i) + vj_uij * deriv_vi_z_env(i)

    jast_elec_Muenv_lapl(i)   = vj_lapl_uij * vi_1b_env(i)                   &
                              + 2.d0 * ( vj_derivx_uij * deriv_vi_x_env(i)   &
                                       + vj_derivy_uij * deriv_vi_y_env(i)   & 
                                       + vj_derivz_uij * deriv_vi_z_env(i) ) &
                              + vj_uij * lapl_vi_env(i)
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, vi_1b_env, (elec_num_8)]

  implicit none
  integer          :: i, iA
  integer          :: ii, phase, b
  double precision :: a, riA, tmp
  double precision :: expo, c

  PROVIDE j1b_type

  if((j1b_type .eq. 2) .or. (j1b_type .eq. 102)) then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = j1b_pen(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp - dexp(-a*riA)
      enddo
      vi_1b_env(i) = tmp
    enddo

  elseif((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = j1b_pen(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp * (1.d0 - dexp(-a*riA*riA))
      enddo
      vi_1b_env(i) = tmp
    enddo

    !do i = 1, elec_num
    !  vi_1b_env(i) = 0.d0
    !  do ii = 1, List_all_comb_b2_size
    !    phase = 0
    !    expo  = 0.d0
    !    !DIR$ LOOP COUNT (100)
    !    do iA = 1, nucl_num
    !      a   = j1b_pen(iA)
    !      b   = List_all_comb_b2(iA,ii)
    !      c   = dble(b) * a
    !      riA = nucl_elec_dist(iA,i)
    !      phase += b
    !      expo  += c * riA * riA
    !    enddo
    !    vi_1b_env(i) += (-1.d0)**dble(phase) * dexp(-expo)
    !  enddo
    !enddo

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = j1b_pen(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp - j1b_pen_coef(iA) * dexp(-a*riA*riA)
      enddo
      vi_1b_env(i) = tmp
    enddo

  elseif((j1b_type .eq. 5) .or. (j1b_type .eq. 105)) then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = j1b_pen(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp - dexp(-a*riA*riA*riA*riA)
      enddo
      vi_1b_env(i) = tmp
    enddo

  else

    print*, 'j1_type = ', j1b_type, 'not implemented yet'
    stop

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, deriv_vi_x_env, (elec_num_8)]
&BEGIN_PROVIDER [double precision, deriv_vi_y_env, (elec_num_8)]
&BEGIN_PROVIDER [double precision, deriv_vi_z_env, (elec_num_8)]
&BEGIN_PROVIDER [double precision,    lapl_vi_env, (elec_num_8)]

  implicit none
  integer          :: i, ii, iA, phase, b
  double precision :: a, riA, dx, dy, dz, r2, r4
  double precision :: tmp, tmpx, tmpy, tmpz, tmpl
  double precision :: expo, coef, coef_x, coef_y, coef_z, c
  double precision :: arg

  PROVIDE j1b_type

  if((j1b_type .eq. 2) .or. (j1b_type .eq. 102)) then

    do i = 1, elec_num
      tmpx = 0.d0
      tmpy = 0.d0
      tmpz = 0.d0
      tmpl = 0.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = j1b_pen(iA)
        ! xi - xA = nucl_elec_dist_vec(1,iA,i)
        dx  = nucl_elec_dist_vec(1,iA,i)
        dy  = nucl_elec_dist_vec(2,iA,i)
        dz  = nucl_elec_dist_vec(3,iA,i)
        riA = nucl_elec_dist(iA,i)
        tmp = a * dexp(-a*riA) / riA
        tmpx = tmpx + tmp * dx
        tmpy = tmpy + tmp * dy
        tmpz = tmpz + tmp * dz
        tmpl = tmpl + tmp * (2.d0 - a*riA)
      enddo
      deriv_vi_x_env(i) = tmpx
      deriv_vi_y_env(i) = tmpy
      deriv_vi_z_env(i) = tmpz
      lapl_vi_env   (i) = tmpl
    enddo

  elseif((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    do i = 1, elec_num
      deriv_vi_x_env(i) = 0.d0
      deriv_vi_y_env(i) = 0.d0
      deriv_vi_z_env(i) = 0.d0
      lapl_vi_env   (i) = 0.d0
      do ii = 1, List_all_comb_b2_size
        phase  = 0
        expo   = 0.d0
        coef   = 0.d0
        coef_x = 0.d0
        coef_y = 0.d0
        coef_z = 0.d0
        !DIR$ LOOP COUNT (100)
        do iA = 1, nucl_num
          a   = j1b_pen(iA)
          b   = List_all_comb_b2(iA,ii)
          c   = dble(b) * a
          riA = nucl_elec_dist(iA,i)
          phase  += b
          coef   += c
          expo   += c * riA * riA
          ! xi - xA = nucl_elec_dist_vec(1,iA,i)
          coef_x += c * nucl_elec_dist_vec(1,iA,i)
          coef_y += c * nucl_elec_dist_vec(2,iA,i)
          coef_z += c * nucl_elec_dist_vec(3,iA,i)
        enddo
        tmp = -2.d0 * (-1.d0)**dble(phase) * dexp(-expo)
        deriv_vi_x_env(i) += tmp * coef_x
        deriv_vi_y_env(i) += tmp * coef_y
        deriv_vi_z_env(i) += tmp * coef_z
        lapl_vi_env   (i) += tmp * (3.d0 * coef - 2.d0 * (coef_x*coef_x + coef_y*coef_y + coef_z*coef_z))
      enddo
    enddo

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    do i = 1, elec_num
      tmpx = 0.d0
      tmpy = 0.d0
      tmpz = 0.d0
      tmpl = 0.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = j1b_pen(iA)
        ! xi - xA = nucl_elec_dist_vec(1,iA,i)
        dx  = nucl_elec_dist_vec(1,iA,i)
        dy  = nucl_elec_dist_vec(2,iA,i)
        dz  = nucl_elec_dist_vec(3,iA,i)
        riA = nucl_elec_dist(iA,i)
        arg = a * riA * riA
        tmp = a * j1b_pen_coef(iA) * dexp(-arg)
        tmpx = tmpx + tmp * dx
        tmpy = tmpy + tmp * dy
        tmpz = tmpz + tmp * dz
        tmpl = tmpl + tmp * (3.d0 - 2.d0 * arg)
      enddo
      deriv_vi_x_env(i) = 2.d0 * tmpx
      deriv_vi_y_env(i) = 2.d0 * tmpy
      deriv_vi_z_env(i) = 2.d0 * tmpz
      lapl_vi_env   (i) = 2.d0 * tmpl
    enddo

  elseif((j1b_type .eq. 5) .or. (j1b_type .eq. 105)) then

    do i = 1, elec_num
      tmpx = 0.d0
      tmpy = 0.d0
      tmpz = 0.d0
      tmpl = 0.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = j1b_pen(iA)
        ! xi - xA = nucl_elec_dist_vec(1,iA,i)
        dx  = nucl_elec_dist_vec(1,iA,i)
        dy  = nucl_elec_dist_vec(2,iA,i)
        dz  = nucl_elec_dist_vec(3,iA,i)
        riA = nucl_elec_dist(iA,i)
        r2  = riA * riA
        r4  = r2  * r2
        tmp = a * r2 * dexp(-a*r4)
        tmpx = tmpx + tmp * dx
        tmpy = tmpy + tmp * dy
        tmpz = tmpz + tmp * dz
        tmpl = tmpl + tmp * (5.d0 - 4.d0*a*r4)
      enddo
      deriv_vi_x_env(i) = 4.d0 * tmpx
      deriv_vi_y_env(i) = 4.d0 * tmpy
      deriv_vi_z_env(i) = 4.d0 * tmpz
      lapl_vi_env   (i) = 4.d0 * tmpl
    enddo

  else

    print*, 'j1_type = ', j1b_type, 'not implemented yet'
    stop

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Muenv_grad_x_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Muenv_grad_y_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Muenv_grad_z_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Muenv_lapl_num  , (elec_num_8)]

  implicit none
  integer          :: i, j
  double precision :: eps, tmp_der, tmp_lap, je_p, je_m, je_0, tmp
  double precision :: r1(3), r2(3)

  eps     = 1d-3
  tmp_der = 0.5d0 /  eps
  tmp_lap = 1.0d0 / (eps * eps)

  do i = 1, elec_num

    r1(1) = elec_coord(i,1)
    r1(2) = elec_coord(i,2)
    r1(3) = elec_coord(i,3)

    jast_elec_Muenv_grad_x_num(i) = 0.d0
    jast_elec_Muenv_grad_y_num(i) = 0.d0
    jast_elec_Muenv_grad_z_num(i) = 0.d0
    jast_elec_Muenv_lapl_num  (i) = 0.d0

    do j = 1, elec_num
      if(j==i) cycle

      r2(1) = elec_coord(j,1)
      r2(2) = elec_coord(j,2)
      r2(3) = elec_coord(j,3)
      call j_elec_Muenv(r1, r2, je_0)

      r1(1) += eps
      call j_elec_Muenv(r1, r2, je_p)
      r1(1) -= 2.d0 * eps
      call j_elec_Muenv(r1, r2, je_m)
      r1(1) += eps
      jast_elec_Muenv_grad_x_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Muenv_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      r1(2) += eps
      call j_elec_Muenv(r1, r2, je_p)
      r1(2) -= 2.d0 * eps
      call j_elec_Muenv(r1, r2, je_m)
      r1(2) += eps
      jast_elec_Muenv_grad_y_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Muenv_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      r1(3) += eps
      call j_elec_Muenv(r1, r2, je_p)
      r1(3) -= 2.d0 * eps
      call j_elec_Muenv(r1, r2, je_m)
      r1(3) += eps
      jast_elec_Muenv_grad_z_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Muenv_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)
    enddo
  enddo

END_PROVIDER

! ---

