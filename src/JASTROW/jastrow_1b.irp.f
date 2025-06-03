! ---

BEGIN_PROVIDER [ double precision, jast_1b_value, (elec_num)  ]

  BEGIN_DOC  
  ! 1-body Jastrow
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: a, b, c, rij, tmp
  double precision :: z, mu, mu_pi, zr, mur

  do i = 1, elec_num

    jast_1b_value(i) = 0.d0

    if( j1b_type .eq. 1 ) then ! add 1body-Gauss Jastrow
                               ! J(i) = - \sum_A [ 1 - exp( -alpha_A r_iA^2 ) ]
      !DIR$ LOOP COUNT (100)
      do j = 1, nucl_num
        a   = j1b_pen(j)
        rij = nucl_elec_dist(j,i)
        tmp = 1.d0 - dexp(-a*rij*rij)
        jast_1b_value(i) -= tmp
      enddo

    elseif( j1b_type .eq. 2 ) then ! add 1body-Gauss Jastrow
                                   ! J(i) = \sum_A [ c_A exp( -alpha_A r_iA^2 ) ]
      !DIR$ LOOP COUNT (100)
      do j = 1, nucl_num
        a   = j1b_pen  (j)
        b   = j1b_coeff(j)
        rij = nucl_elec_dist(j,i)
        tmp = b * dexp(-a*rij*rij)
        jast_1b_value(i) += tmp
      enddo

!    elseif( j1b_type .eq. 3 ) then ! add 1body-Slater Jastrow
!                                   ! J(i) = - \sum_A c_A exp( - alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1bslat_expo(j)
!        c   = jast_1bslat_coef(j)
!        rij = nucl_elec_dist(j,i)
!        tmp = c * dexp( - a * rij )
!        jast_1b_value(i) -= tmp
!      enddo
!
!    elseif( j1b_type .eq. 4 ) then ! add 1body-Tanh Jastrow
!                                   ! J(i) = - \sum_A tanh(alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1btanh_pen(j)
!        rij = nucl_elec_dist(j,i)
!        tmp = dtanh(a*rij)
!        jast_1b_value(i) -= tmp
!      enddo
!
!    elseif( j1b_type .eq. 5 ) then ! add 1body-Simple Jastrow
!                                   ! J(i) = - \sum_A [ (alpha_A r_iA) / (1 + alpha_A r_iA) ]^2
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_pen(j)
!        rij = a * nucl_elec_dist(j,i)
!        tmp = rij / (1.d0 + rij)
!        jast_1b_value(i) -= tmp*tmp
!      enddo
!
!    elseif( j1b_type .eq. 6 ) then ! add 1body-RSDFT  Jastrow
!                                   ! J(i) = - \sum_A [ -z_A r_iA erfc(mu*r_iA) + z_A exp(-(mu*r_iA)^2)/(mu*sqt_pi) ]
!      mu    = mu_erf
!      mu_pi = 1.d0 / ( dsqpi * mu )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        rij = nucl_elec_dist(j,i)
!        z   = nucl_charge(j) 
!        zr  = z  * rij
!        mur = mu * rij
!        tmp = - zr * ( 1.d0 - derf(mur) ) + z * mu_pi * dexp(-mur*mur)
!        jast_1b_value(i) -= tmp
!      enddo
!
!    elseif( j1b_type .eq. 7 ) then ! add 1body-erf Jastrow
!                                   ! J(i) = - \sum_A erf( alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1berf_pen(j)
!        rij = nucl_elec_dist(j,i)
!        tmp = derf(a*rij)
!        jast_1b_value(i) -= tmp
!      enddo

    endif
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, jast_1b_grad_x, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_1b_grad_y, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_1b_grad_z, (elec_num) ]

  BEGIN_DOC  
  ! Gradient of the Jastrow
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: a, b, c, rij, tmp
  double precision :: z, mu, mur

  do i = 1, elec_num

    jast_1b_grad_x(i) = 0.d0
    jast_1b_grad_y(i) = 0.d0
    jast_1b_grad_z(i) = 0.d0


    if( j1b_type .eq. 1 ) then ! add 1body-Gauss Jastrow
                               ! J(i) = - \sum_A [ 1 - exp( -alpha_A r_iA^2 ) ]
      !DIR$ LOOP COUNT (100)
      do j = 1, nucl_num
        a   = j1b_pen(j)
        rij = nucl_elec_dist(j,i)
        tmp = 2.d0 * a * dexp(-a*rij*rij)
        jast_1b_grad_x(i) -= nucl_elec_dist_vec(1,j,i) * tmp
        jast_1b_grad_y(i) -= nucl_elec_dist_vec(2,j,i) * tmp
        jast_1b_grad_z(i) -= nucl_elec_dist_vec(3,j,i) * tmp
      enddo

    elseif( j1b_type .eq. 2 ) then ! add 1body-Gauss Jastrow
                                   ! J(i) = \sum_A [ c_A exp( -alpha_A r_iA^2 ) ]
      !DIR$ LOOP COUNT (100)
      do j = 1, nucl_num
        a   = j1b_pen  (j)
        b   = j1b_coeff(j)
        rij = nucl_elec_dist(j,i)
        tmp = 2.d0 * a * b * dexp(-a*rij*rij)
        jast_1b_grad_x(i) -= nucl_elec_dist_vec(1,j,i) * tmp
        jast_1b_grad_y(i) -= nucl_elec_dist_vec(2,j,i) * tmp
        jast_1b_grad_z(i) -= nucl_elec_dist_vec(3,j,i) * tmp
      enddo

!    elseif( j1b_type .eq. 3 ) then ! add 1body-Slater Jastrow
!                                   ! J(i) = - \sum_A c_A exp( - alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1bslat_expo(j)
!        c   = jast_1bslat_coef(j)
!        rij = nucl_elec_dist(j,i)
!        tmp = c * a * dexp( - a * rij ) / rij
!        jast_1b_grad_x(i) -= nucl_elec_dist_vec(1,j,i) * tmp
!        jast_1b_grad_y(i) -= nucl_elec_dist_vec(2,j,i) * tmp
!        jast_1b_grad_z(i) -= nucl_elec_dist_vec(3,j,i) * tmp
!      enddo
!
!    elseif( j1b_type .eq. 4 ) then ! add 1body-Tanh Jastrow
!                                   ! J(i) = - \sum_A tanh(alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1btanh_pen(j)
!        rij = nucl_elec_dist(j,i)
!        c   = dtanh(a*rij)
!        tmp = a * ( 1.d0 - c*c ) / rij
!        jast_1b_grad_x(i) -= nucl_elec_dist_vec(1,j,i) * tmp
!        jast_1b_grad_y(i) -= nucl_elec_dist_vec(2,j,i) * tmp
!        jast_1b_grad_z(i) -= nucl_elec_dist_vec(3,j,i) * tmp
!      enddo
!
!    elseif( j1b_type .eq. 5 ) then ! add 1body-Simple Jastrow
!                                   ! J(i) = - \sum_A [ (alpha_A r_iA) / (1 + alpha_A r_iA) ]^2
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_pen(j)
!        rij = a * nucl_elec_dist(j,i)
!        tmp = (a+a)*a / (1.d0+rij*(3.d0+rij*(3.d0+rij)))
!        jast_1b_grad_x(i) -= nucl_elec_dist_vec(1,j,i) * tmp
!        jast_1b_grad_y(i) -= nucl_elec_dist_vec(2,j,i) * tmp
!        jast_1b_grad_z(i) -= nucl_elec_dist_vec(3,j,i) * tmp
!      enddo
!
!    elseif( j1b_type .eq. 6 ) then ! add 1body-RSDFT  Jastrow
!                                   ! J(i) = - \sum_A [ -z_A r_iA erfc(mu*r_iA) + z_A exp(-(mu*r_iA)^2)/(mu*sqt_pi) ]
!      mu = mu_erf
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        rij = nucl_elec_dist(j,i)
!        z   = nucl_charge(j) 
!        mur = mu * rij
!        tmp = -z * ( 1.d0 - derf(mur) ) / rij 
!        jast_1b_grad_x(i) -= nucl_elec_dist_vec(1,j,i) * tmp
!        jast_1b_grad_y(i) -= nucl_elec_dist_vec(2,j,i) * tmp
!        jast_1b_grad_z(i) -= nucl_elec_dist_vec(3,j,i) * tmp
!      enddo
!
!    elseif( j1b_type .eq. 7 ) then ! add 1body-erf Jastrow
!                                   ! J(i) = - \sum_A erf( alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1berf_pen(j)
!        rij = nucl_elec_dist(j,i)
!        c   = a * rij
!        tmp = 2.d0 * a * dexp(-c*c) / (dsqpi * rij)
!        jast_1b_grad_x(i) -= nucl_elec_dist_vec(1,j,i) * tmp
!        jast_1b_grad_y(i) -= nucl_elec_dist_vec(2,j,i) * tmp
!        jast_1b_grad_z(i) -= nucl_elec_dist_vec(3,j,i) * tmp
!      enddo

    endif
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_1b_lapl, (elec_num) ]

  BEGIN_DOC  
  ! Laplacian of the Jastrow factor
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: a, b, c, rij, tmp
  double precision :: mu, mu_pi, mur, z

  do i = 1, elec_num

    jast_1b_lapl(i) = 0.d0

    if( j1b_type .eq. 1 ) then ! add 1body-Gauss Jastrow
                               ! J(i) = - \sum_A [ 1 - exp( -alpha_A r_iA^2 ) ]
      !DIR$ LOOP COUNT (100)
      do j = 1, nucl_num
        a   = j1b_pen(j)
        rij = nucl_elec_dist(j,i)
        c   = a * rij * rij
        tmp = 2.d0 * a * dexp(-c) * (3.d0-2.d0*c)
        jast_1b_lapl(i) -= tmp
      enddo

    elseif( j1b_type .eq. 2 ) then ! add 1body-Gauss Jastrow
                                   ! J(i) = \sum_A [ c_A exp( -alpha_A r_iA^2 ) ]
      !DIR$ LOOP COUNT (100)
      do j = 1, nucl_num
        a   = j1b_pen  (j)
        b   = j1b_coeff(j)
        rij = nucl_elec_dist(j,i)
        c   = a * rij * rij
        tmp = 2.d0 * a * b * dexp(-c) * (3.d0-2.d0*c)
        jast_1b_lapl(i) -= tmp
      enddo


!    elseif( j1b_type .eq. 3 ) then ! add 1body-Slater Jastrow
!                                       ! J(i) = - \sum_A c_A exp( - alpha_A r_iA )
!    !DIR$ LOOP COUNT (100)
!    do j = 1, nucl_num
!      a   = jast_1bslat_expo(j)
!      c   = jast_1bslat_coef(j)
!      rij = nucl_elec_dist(j,i)
!      tmp = c * a * dexp(-a*rij) * ( 2.d0/rij - a )
!      jast_1b_lapl(i) -= tmp
!    enddo
!
!    elseif( j1b_type .eq. 4 ) then ! add 1body-Tanh Jastrow
!                                       ! J(i) = - \sum_A tanh(alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1btanh_pen(j)
!        rij = nucl_elec_dist(j,i)
!        c   = dtanh(a*rij)
!        tmp = 2.d0 * a * ( 1.d0 - c*c ) * ( 1.d0/rij - a*c )
!        jast_1b_lapl(i) -= tmp
!      enddo
!
!    elseif( j1b_type .eq. 5 ) then ! add 1body-Simple Jastrow
!                                       ! J(i) = - \sum_A [ (alpha_A r_iA) / (1 + alpha_A r_iA) ]^2
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_pen(j)
!        rij = a * nucl_elec_dist(j,i)
!        tmp = 6.d0*a*a / (1.d0+rij*(4.d0+rij*(6.d0+rij*(4.d0+rij))))
!        jast_1b_lapl(i) -= tmp
!      enddo
!
!    elseif( j1b_type .eq. 6 ) then ! add 1body-RSDFT  Jastrow
!                                       ! J(i) = - \sum_A [ -z_A r_iA erfc(mu*r_iA) + z_A exp(-(mu*r_iA)^2)/(mu*sqt_pi) ]
!      mu    = mu_erf
!      mu_pi = mu / dsqpi
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        rij = nucl_elec_dist(j,i)
!        z   = nucl_charge(j) 
!        mur = mu * rij
!        tmp = -2.d0*z*(1.d0-derf(mur))/rij + 2.d0*z*mu_pi*dexp(-mur*mur)
!        jast_1b_lapl(i) -= tmp
!      enddo
!
!    elseif( j1b_type .eq. 7 ) then ! add 1body-erf Jastrow
!                                       ! J(i) = - \sum_A erf( alpha_A r_iA )
!      !DIR$ LOOP COUNT (100)
!      do j = 1, nucl_num
!        a   = jast_1berf_pen(j)
!        rij = nucl_elec_dist(j,i)
!        c   = a * rij
!        tmp = 4.d0 * dexp(-c*c) * (a/rij-a*a*a*rij) / dsqpi
!        jast_1b_lapl(i) -= tmp
!      enddo

    endif
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_1b_grad_sq, (elec_num) ]

  BEGIN_DOC  
  ! square of the gradient of the 1-body Jastrow
  END_DOC

  implicit none
  integer :: i

  do i = 1, elec_num
    jast_1b_grad_sq(i) = jast_1b_grad_x(i) * jast_1b_grad_x(i) &
                       + jast_1b_grad_y(i) * jast_1b_grad_y(i) &
                       + jast_1b_grad_z(i) * jast_1b_grad_z(i) 
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_J1b_grad ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_1b_grad_sq(i)
  enddo
  deltaE_J1b_grad = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_J1b_lapl ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_1b_lapl(i)
  enddo
  deltaE_J1b_lapl = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_J1b_nonh ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= psidet_right_grad_lapl(1,i) * jast_1b_grad_x(i) &
         + psidet_right_grad_lapl(2,i) * jast_1b_grad_y(i) &
         + psidet_right_grad_lapl(3,i) * jast_1b_grad_z(i)
  enddo
  deltaE_J1b_nonh = tmp * psidet_right_inv

END_PROVIDER

! ---

