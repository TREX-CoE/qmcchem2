! Mu Jastrow
! --------------

 BEGIN_PROVIDER [ double precision, jast_Mu_value ]
&BEGIN_PROVIDER [ double precision, jast_Mu_value_inv ]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Mu_value(i)
  enddo

  jast_Mu_value     = dexp(argexpo)
  jast_Mu_value_inv = 1.d0 / jast_Mu_value

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Mu_value, (elec_num_8)  ]

  BEGIN_DOC  
  ! J(i) = 0.50 \sum_{j!=i} 0.50 [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ]
  !      = 0.25 \sum_{j!=i}      [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ]
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: rij, tmp, mu, mu_sq, mu_pi

  mu    = mu_erf
  mu_sq = mu * mu
  mu_pi = 1.d0 / ( dsqpi * mu )

  do i = 1, elec_num
    jast_elec_Mu_value(i) = 0.d0
    !DIR$ LOOP COUNT(100)
    do j = 1, elec_num 
      if(j==i) cycle
      rij = elec_dist(j,i)
      tmp = 0.25d0 * ( rij*(1.d0-derf(mu*rij)) - mu_pi*dexp(-mu_sq*rij*rij) )
      jast_elec_Mu_value(i) += tmp 
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision , jast_elec_Mu_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Mu_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Mu_grad_z, (elec_num_8) ]

  BEGIN_DOC  
  ! Gradient of the Jastrow factor
  ! eq (A1) in M. Giner, JCP 2021
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: rij, mu, tmp, xx, yy, zz

  mu = mu_erf

  do i = 1, elec_num
    jast_elec_Mu_grad_x(i) = 0.d0
    jast_elec_Mu_grad_y(i) = 0.d0
    jast_elec_Mu_grad_z(i) = 0.d0
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num
      if(i==j) cycle
      rij = elec_dist(j,i)
      tmp = 0.5d0 * (derf(mu * rij) - 1.d0) 

      xx  = elec_dist_inv(j,i) * elec_dist_vec_x(j,i)
      yy  = elec_dist_inv(j,i) * elec_dist_vec_y(j,i)
      zz  = elec_dist_inv(j,i) * elec_dist_vec_z(j,i)

      jast_elec_Mu_grad_x(i) += tmp * xx 
      jast_elec_Mu_grad_y(i) += tmp * yy 
      jast_elec_Mu_grad_z(i) += tmp * zz 
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Mu_lapl, (elec_num_8) ]

  BEGIN_DOC  
  ! Laplacian of the Jastrow factor
  !      [1 - erf(mu r)]/r - mu exp[-(mu r)^2]/sqrt(pi)
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: rij, mu, mu_pi

  mu    = mu_erf
  mu_pi = mu / dsqpi

  do i = 1, elec_num
    jast_elec_Mu_lapl(i) = 0.d0
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num
      if(j==i) cycle
      rij = mu * elec_dist(j,i)
      jast_elec_Mu_lapl(i) += (1.d0 - derf(rij)) * elec_dist_inv(j,i) - mu_pi * dexp(-rij*rij)
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Mu_lapl_reg, (elec_num_8) ]

  BEGIN_DOC  
  ! Regularized laplacian of the Jastrow:
  !      [1 - erf(mu r)]/r - mu exp[-(mu r)^2]/sqrt(pi) - 1/r 
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: rij, mu, mu_pi, tmp

  mu    = mu_erf
  mu_pi = mu / dsqpi
  do i = 1, elec_num
    jast_elec_Mu_lapl_reg(i) = 0.d0
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num
      if(j==i) cycle
      rij = mu * elec_dist(j,i)
      tmp = derf(rij) * elec_dist_inv(j,i) + mu_pi * dexp(-rij*rij)
      jast_elec_Mu_lapl_reg(i) -= tmp 
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, grad_j_mu_x,(elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_y,(elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_z,(elec_num, elec_num)]

  BEGIN_DOC  
  ! useful for three_body_mu calculation 
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: rij, mu, tmp, xx, yy, zz

  mu = mu_erf
  grad_j_mu_x = 0.d0
  grad_j_mu_y = 0.d0
  grad_j_mu_z = 0.d0

  do j = 1, elec_num
    do i = 1, elec_num
      if(i==j) cycle
      rij = elec_dist(i,j)
      tmp = 0.5d0 * ( 1.d0 - derf(mu * rij) ) 

      xx = (elec_coord_transp(1,i) - elec_coord_transp(1,j)) * elec_dist_inv(i,j)
      yy = (elec_coord_transp(2,i) - elec_coord_transp(2,j)) * elec_dist_inv(i,j)
      zz = (elec_coord_transp(3,i) - elec_coord_transp(3,j)) * elec_dist_inv(i,j)

      grad_j_mu_x(i,j) = xx * tmp 
      grad_j_mu_y(i,j) = yy * tmp
      grad_j_mu_z(i,j) = zz * tmp
    enddo
  enddo

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, grad_j_mu_coul_x, (elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_coul_y, (elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_coul_z, (elec_num, elec_num)]

  implicit none
  integer          :: i, j
  double precision :: rij, mu, scal

  mu = mu_erf
  grad_j_mu_coul_x = 0.d0
  grad_j_mu_coul_y = 0.d0
  grad_j_mu_coul_z = 0.d0
  do j = 1, elec_num
    do i = 1, elec_num
      if(i==j) cycle
      rij  = elec_dist(i,j)
      scal = 0.5d0 * elec_dist_inv(i,j) 
      grad_j_mu_coul_x(i,j) = (elec_coord_transp(1,i) - elec_coord_transp(1,j)) * scal
      grad_j_mu_coul_y(i,j) = (elec_coord_transp(2,i) - elec_coord_transp(2,j)) * scal
      grad_j_mu_coul_z(i,j) = (elec_coord_transp(3,i) - elec_coord_transp(3,j)) * scal
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu_lapl ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_lapl(i)
  enddo
  deltaE_Jmu_lapl = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu_lapl_reg ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_lapl_reg(i)
  enddo
  deltaE_Jmu_lapl_reg = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu_nonh ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= ( psidet_right_grad_lapl(1,i) * jast_elec_Mu_grad_x(i) &
           + psidet_right_grad_lapl(2,i) * jast_elec_Mu_grad_y(i) &
           + psidet_right_grad_lapl(3,i) * jast_elec_Mu_grad_z(i) ) * psidet_right_inv
  enddo
  deltaE_Jmu_nonh = tmp 

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu_grad ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_grad_x(i) * jast_elec_Mu_grad_x(i) &
         + jast_elec_Mu_grad_y(i) * jast_elec_Mu_grad_y(i) &
         + jast_elec_Mu_grad_z(i) * jast_elec_Mu_grad_z(i)
  enddo
  deltaE_Jmu_grad = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu_grad_2b ]

  implicit none
  integer          :: i, j
  double precision :: mu, tmp, rij, a

  mu  = mu_erf
  tmp = 0.d0
  do i = 1, elec_num
    !DIR$ LOOP COUNT(100)
    do j = 1, elec_num
      if(i==j) cycle
      rij = mu * elec_dist(j,i)
      a   = 1.d0 - derf(rij)
      tmp = tmp - a * a
    enddo
  enddo

  ! x 0.50 for double couting
  ! x 0.25 formula
  deltaE_Jmu_grad_2b = 0.125d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu_grad_3b ]

  implicit none
  integer :: i, j, k

  deltaE_Jmu_grad_3b = 0.d0
  do i = 1, elec_num
    do j = i+1, elec_num
      do k = j+1, elec_num
        deltaE_Jmu_grad_3b -= grad_j_mu_x(i,j) * grad_j_mu_x(i,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_y(i,j) * grad_j_mu_y(i,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_z(i,j) * grad_j_mu_z(i,k)

        deltaE_Jmu_grad_3b -= grad_j_mu_x(j,i) * grad_j_mu_x(j,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_y(j,i) * grad_j_mu_y(j,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_z(j,i) * grad_j_mu_z(j,k)

        deltaE_Jmu_grad_3b -= grad_j_mu_x(k,i) * grad_j_mu_x(k,j)
        deltaE_Jmu_grad_3b -= grad_j_mu_y(k,i) * grad_j_mu_y(k,j)
        deltaE_Jmu_grad_3b -= grad_j_mu_z(k,i) * grad_j_mu_z(k,j)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

