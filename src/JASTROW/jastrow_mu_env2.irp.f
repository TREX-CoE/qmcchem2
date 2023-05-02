! Mu Jastrow x envelop_type 2
! ---------------------------

! ---

subroutine v1b_env2(i, j1b)

  implicit none
  integer,          intent(in)  :: i
  double precision, intent(out) :: j1b
  integer                       :: iA
  double precision              :: a, riA

  j1b = 1.d0
  !DIR$ LOOP COUNT (100)
  do iA = 1, nucl_num
    a   = j1b_pen(iA)
    riA = nucl_elec_dist(iA,i)
    j1b = j1b - dexp(-a*riA)
  enddo

  return
end subroutine v1b_env2

! ---

BEGIN_PROVIDER [double precision, vi_1b_env2, (elec_num_8)]

  implicit none
  integer          :: i, iA
  double precision :: a, riA, tmp

  do i = 1, elec_num

    tmp = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = j1b_pen(iA)
      riA = nucl_elec_dist(iA,i)

      tmp = tmp - dexp(-a*riA)
    enddo

    vi_1b_env2(i) = tmp
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, deriv_vi_x_env2, (elec_num_8)]
&BEGIN_PROVIDER [double precision, deriv_vi_y_env2, (elec_num_8)]
&BEGIN_PROVIDER [double precision, deriv_vi_z_env2, (elec_num_8)]
&BEGIN_PROVIDER [double precision,    lapl_vi_env2, (elec_num_8)]

  implicit none
  integer          :: i, ii, iA, phase, b
  double precision :: a, riA, dx, dy, dz
  double precision :: tmp, tmpx, tmpy, tmpz, tmpl

  do i = 1, elec_num

    tmpx = 0.d0
    tmpy = 0.d0
    tmpz = 0.d0
    tmpl = 0.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = j1b_pen(iA)

      ! xi - xA = nucl_elec_dist_vec(1,iA,i)
      dx = nucl_elec_dist_vec(1,iA,i)
      dy = nucl_elec_dist_vec(2,iA,i)
      dz = nucl_elec_dist_vec(3,iA,i)

      riA = dsqrt(dx*dx + dy*dy + dz*dz)
      tmp = a * dexp(-a*riA) / riA

      tmpx = tmpx + tmp * dx
      tmpy = tmpy + tmp * dy
      tmpz = tmpz + tmp * dz
      tmpl = tmpl + tmp * (2.d0 - a*riA)
    enddo

    deriv_vi_x_env2(i) = tmpx
    deriv_vi_y_env2(i) = tmpy
    deriv_vi_z_env2(i) = tmpz
    lapl_vi_env2   (i) = tmpl
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_Mu_env2_value ]
&BEGIN_PROVIDER [double precision, jast_Mu_env2_value_inv ]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Mu_env2_value(i)
  enddo

  jast_Mu_env2_value     = dexp(argexpo)
  jast_Mu_env2_value_inv = 1.d0 / jast_Mu_env2_value

END_PROVIDER

! ---

subroutine j_elec_mu_env2(i, j, je)

  BEGIN_DOC  
  !
  ! J(i,j) = 0.5 [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ] 
  !        x [1 - \sum_A exp(-a_A riA)]
  !        x [1 - \sum_A exp(-a_A rjA)]
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer,          intent(in)  :: i, j
  double precision, intent(out) :: je
  double precision              :: mu_rij
  double precision              :: rij, u_ij, vi, vj

  rij    = elec_dist(j,i)
  mu_rij = mu_erf * rij
  u_ij   = 0.5d0 * (rij * (1.d0 - derf(mu_rij)) - dexp(-mu_rij*mu_rij)/(dsqpi*mu_erf))

  call v1b_env2(i, vi)
  call v1b_env2(j, vj)

  je = u_ij * vi * vj

  return
end subroutine j_elec_mu_env2

! ---

BEGIN_PROVIDER [double precision, jast_elec_Mu_env2_value, (elec_num_8)]

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

  mu_pi = 1.d0 / ( dsqpi * mu_erf )

  do i = 1, elec_num

    tmp_ij = 0.d0
    !DIR$ LOOP COUNT(100)
    do j = 1, elec_num

      if(j==i) cycle

      rij    = elec_dist(j,i)
      mu_rij = mu_erf * rij

      tmp_ij += (rij * (1.d0 - derf(mu_rij)) - mu_pi * dexp(-mu_rij*mu_rij)) * vi_1b_env2(j)
    enddo

    jast_elec_Mu_env2_value(i) = 0.25d0 * tmp_ij * vi_1b_env2(i)
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, jast_elec_Mu_env2_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env2_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env2_grad_z, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env2_lapl  , (elec_num_8) ]

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
      tmp3_ij = -0.5d0 * tmp2_ij * vi_1b_env2(j)

      vj_uij        += 0.5d0 * (rij * tmp1_ij - mu_sqrtpi_inv * tmp0_ij) * vi_1b_env2(j)
      vj_derivx_uij += tmp3_ij * elec_dist_vec_x(j,i)
      vj_derivy_uij += tmp3_ij * elec_dist_vec_y(j,i)
      vj_derivz_uij += tmp3_ij * elec_dist_vec_z(j,i)
      vj_lapl_uij   += (tmp2_ij - mu_div_sqrtpi * tmp0_ij) * vi_1b_env2(j)
    enddo

    jast_elec_Mu_env2_grad_x(i) = vj_derivx_uij * vi_1b_env2(i) + vj_uij * deriv_vi_x_env2(i)
    jast_elec_Mu_env2_grad_y(i) = vj_derivy_uij * vi_1b_env2(i) + vj_uij * deriv_vi_y_env2(i)
    jast_elec_Mu_env2_grad_z(i) = vj_derivz_uij * vi_1b_env2(i) + vj_uij * deriv_vi_z_env2(i)

    jast_elec_Mu_env2_lapl(i)   = vj_lapl_uij * vi_1b_env2(i)                    &
                                 + 2.d0 * ( vj_derivx_uij * deriv_vi_x_env2(i)   &
                                          + vj_derivy_uij * deriv_vi_y_env2(i)   & 
                                          + vj_derivz_uij * deriv_vi_z_env2(i) ) &
                                 + vj_uij * lapl_vi_env2(i)
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Mu_env2_grad_x_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_env2_grad_y_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_env2_grad_z_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_env2_lapl_num  , (elec_num_8)]

  implicit none
  integer          :: i, j
  double precision :: eps, tmp_der, tmp_lap, je_p, je_m, je_0, tmp

  eps     = 1d-3
  tmp_der = 0.5d0 /  eps
  tmp_lap = 1.0d0 / (eps * eps)

  do i = 1, elec_num

    ! --- --- ---
    ! d / dx

    tmp = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      elec_coord(i,1) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_p)
      
      elec_coord(i,1) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_m)
      
      elec_coord(i,1) += eps
      TOUCH elec_coord

      tmp += tmp_der * (je_p - je_m)
    enddo

    jast_elec_Mu_env2_grad_x_num(i) = tmp
    !
    ! --- --- ---

    ! --- --- ---
    ! d / dy

    tmp = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      elec_coord(i,2) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_p)
      
      elec_coord(i,2) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_m)
      
      elec_coord(i,2) += eps
      TOUCH elec_coord

      tmp += tmp_der * (je_p - je_m)
    enddo

    jast_elec_Mu_env2_grad_y_num(i) = tmp
    !
    ! --- --- ---

    ! --- --- ---
    ! d / dz

    tmp = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      elec_coord(i,3) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_p)

      elec_coord(i,3) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_m)

      elec_coord(i,3) += eps
      TOUCH elec_coord

      tmp += tmp_der * (je_p - je_m)
    enddo

    jast_elec_Mu_env2_grad_z_num(i) = tmp
    !
    ! --- --- ---

    ! --- --- ---
    ! d^2 / dx^2 + d^2 / dy^2 + d^2 / dz^2 

    tmp = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      elec_coord(i,1) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_p)
      
      elec_coord(i,1) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_m)
      
      elec_coord(i,1) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_0)

      tmp += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      ! ---

      elec_coord(i,2) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_p)
      
      elec_coord(i,2) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_m)
      
      elec_coord(i,2) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_0)
      
      tmp += tmp_lap * (je_p - 2.d0 * je_0 + je_m)
      
      ! ---
      
      elec_coord(i,3) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_p)
      
      elec_coord(i,3) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_m)
      
      elec_coord(i,3) += eps
      TOUCH elec_coord
      call j_elec_mu_env2(i, j, je_0)

      tmp += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

    enddo

    jast_elec_Mu_env2_lapl_num(i) = tmp
    !
    ! --- --- ---

  enddo

END_PROVIDER

