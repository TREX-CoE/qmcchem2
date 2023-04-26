! Mu Jastrow x envelop_type 3
! ---------------------------

! ---

BEGIN_PROVIDER [ integer, List_all_comb_b2_size]

  implicit none

  List_all_comb_b2_size = 2**nucl_num

END_PROVIDER

! ---

BEGIN_PROVIDER [ integer, List_all_comb_b2, (nucl_num, List_all_comb_b2_size)]

  implicit none
  integer :: i, j

  if(nucl_num .gt. 32) then
    print *, ' nucl_num = ', nucl_num, '> 32'
    stop
  endif

  List_all_comb_b2 = 0

  do i = 0, List_all_comb_b2_size-1
    do j = 0, nucl_num-1
      if (btest(i,j)) then
        List_all_comb_b2(j+1,i+1) = 1
      endif
    enddo
  enddo

END_PROVIDER

! ---

subroutine v1b_env3(i, j1b)

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
    j1b = j1b * (1.d0 - dexp(-a*riA*riA))
  enddo

  return
end subroutine v1b_env3

! ---

BEGIN_PROVIDER [ double precision , vi_1b, (elec_num_8) ]

  implicit none
  integer          :: i, iA
  double precision :: a, riA, tmp

  do i = 1, elec_num

    tmp = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = j1b_pen(iA)
      riA = nucl_elec_dist(iA,i)

      tmp = tmp * (1.d0 - dexp(-a*riA*riA))
      !print*, a, riA, 1.d0 - dexp(-a*riA*riA), tmp
    enddo

    vi_1b(i) = tmp
  enddo

END_PROVIDER

!BEGIN_PROVIDER [ double precision , vi_1b, (elec_num_8) ]
!
!  implicit none
!  integer          :: i, ii, iA, phase, b
!  double precision :: expo, a, c, riA
!
!  do i = 1, elec_num
!
!    vi_1b(i) = 0.d0
!    do ii = 1, List_all_comb_b2_size
!  
!      phase = 0
!      expo  = 0.d0
!      !DIR$ LOOP COUNT (100)
!      do iA = 1, nucl_num
!        a   = j1b_pen(iA)
!        b   = List_all_comb_b2(iA,ii)
!        c   = dble(b) * a
!        riA = nucl_elec_dist(iA,i)
!
!        phase += b
!        expo  += c * riA * riA
!      enddo
!
!      vi_1b(i) += (-1.d0)**dble(phase) * dexp(-expo)
!    enddo
!  enddo
!
!END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision , deriv_vi_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , deriv_vi_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , deriv_vi_z, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision ,    lapl_vi, (elec_num_8) ]

  implicit none
  integer          :: i, ii, iA, phase, b
  double precision :: expo, coef, coef_x, coef_y, coef_z, a, c, riA
  double precision :: tmp

  do i = 1, elec_num

    deriv_vi_x(i) = 0.d0
    deriv_vi_y(i) = 0.d0
    deriv_vi_z(i) = 0.d0
    lapl_vi   (i) = 0.d0
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

      deriv_vi_x(i) += tmp * coef_x 
      deriv_vi_y(i) += tmp * coef_y 
      deriv_vi_z(i) += tmp * coef_z 
      lapl_vi   (i) += tmp * (3.d0 * coef - 2.d0 * (coef_x*coef_x + coef_y*coef_y + coef_z*coef_z))
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, jast_Mu_env3_value ]
&BEGIN_PROVIDER [ double precision, jast_Mu_env3_value_inv ]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Mu_env3_value(i)
  enddo

  jast_Mu_env3_value     = dexp(argexpo)
  jast_Mu_env3_value_inv = 1.d0 / jast_Mu_env3_value

END_PROVIDER

! ---

subroutine j_elec_mu_env3(i, j, je)

  BEGIN_DOC  
  !
  ! J(i,j) = 0.5 [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ] 
  !        x \prod_A [(1 - exp(-a_A riA**2))]
  !        x \prod_A [(1 - exp(-a_A rjA**2))]
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

  call v1b_env3(i, vi)
  call v1b_env3(j, vj)

  je = u_ij * vi * vj

  return
end subroutine j_elec_mu_env3

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Mu_env3_value, (elec_num_8)  ]

  BEGIN_DOC  
  !
  ! J(i) = 0.5  \sum_{j!=i} 0.5 [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ]
  !      x \prod_A [(1 - exp(-a_A riA**2))] [(1 - exp(-a_A rjA**2))]
  !
  !      = 0.25 \sum_{j!=i} [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ]
  !      x \prod_A [(1 - exp(-a_A riA**2))] [(1 - exp(-a_A rjA**2))]
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

      tmp_ij += (rij * (1.d0 - derf(mu_rij)) - mu_pi * dexp(-mu_rij*mu_rij)) * vi_1b(j)
    enddo

    jast_elec_Mu_env3_value(i) = 0.25d0 * tmp_ij * vi_1b(i)
    !print *, i, tmp_ij, vi_1b(i), jast_elec_Mu_env3_value(i)
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_grad_z, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_lapl  , (elec_num_8) ]

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
      tmp3_ij = -0.5d0 * tmp2_ij * vi_1b(j)

      vj_uij        += 0.5d0 * (rij * tmp1_ij - mu_sqrtpi_inv * tmp0_ij) * vi_1b(j)
      vj_derivx_uij += tmp3_ij * elec_dist_vec_x(j,i)
      vj_derivy_uij += tmp3_ij * elec_dist_vec_y(j,i)
      vj_derivz_uij += tmp3_ij * elec_dist_vec_z(j,i)
      vj_lapl_uij   += (tmp2_ij - mu_div_sqrtpi * tmp0_ij) * vi_1b(j)
    enddo

    jast_elec_Mu_env3_grad_x(i) = vj_derivx_uij * vi_1b(i) + vj_uij * deriv_vi_x(i)
    jast_elec_Mu_env3_grad_y(i) = vj_derivy_uij * vi_1b(i) + vj_uij * deriv_vi_y(i)
    jast_elec_Mu_env3_grad_z(i) = vj_derivz_uij * vi_1b(i) + vj_uij * deriv_vi_z(i)

    jast_elec_Mu_env3_lapl(i)   = vj_lapl_uij * vi_1b(i)                   &
                                 + 2.d0 * ( vj_derivx_uij * deriv_vi_x(i)   &
                                          + vj_derivy_uij * deriv_vi_y(i)   & 
                                          + vj_derivz_uij * deriv_vi_z(i) ) &
                                 + vj_uij * lapl_vi(i)
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_grad_x_num, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_grad_y_num, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_grad_z_num, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Mu_env3_lapl_num  , (elec_num_8) ]

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
      call j_elec_mu_env3(i, j, je_p)
      
      elec_coord(i,1) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_m)
      
      elec_coord(i,1) += eps
      TOUCH elec_coord

      tmp += tmp_der * (je_p - je_m)
    enddo

    jast_elec_Mu_env3_grad_x_num(i) = tmp
    !
    ! --- --- ---

    ! --- --- ---
    ! d / dy

    tmp = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      elec_coord(i,2) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_p)
      
      elec_coord(i,2) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_m)
      
      elec_coord(i,2) += eps
      TOUCH elec_coord

      tmp += tmp_der * (je_p - je_m)
    enddo

    jast_elec_Mu_env3_grad_y_num(i) = tmp
    !
    ! --- --- ---

    ! --- --- ---
    ! d / dz

    tmp = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      elec_coord(i,3) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_p)

      elec_coord(i,3) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_m)

      elec_coord(i,3) += eps
      TOUCH elec_coord

      tmp += tmp_der * (je_p - je_m)
    enddo

    jast_elec_Mu_env3_grad_z_num(i) = tmp
    !
    ! --- --- ---

    ! --- --- ---
    ! d^2 / dx^2 + d^2 / dy^2 + d^2 / dz^2 

    tmp = 0.d0
    do j = 1, elec_num
      if(j==i) cycle

      elec_coord(i,1) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_p)
      
      elec_coord(i,1) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_m)
      
      elec_coord(i,1) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_0)

      tmp += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      ! ---

      elec_coord(i,2) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_p)
      
      elec_coord(i,2) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_m)
      
      elec_coord(i,2) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_0)
      
      tmp += tmp_lap * (je_p - 2.d0 * je_0 + je_m)
      
      ! ---
      
      elec_coord(i,3) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_p)
      
      elec_coord(i,3) -= 2.d0 * eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_m)
      
      elec_coord(i,3) += eps
      TOUCH elec_coord
      call j_elec_mu_env3(i, j, je_0)

      tmp += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

    enddo

    jast_elec_Mu_env3_lapl_num(i) = tmp
    !
    ! --- --- ---

  enddo

END_PROVIDER

