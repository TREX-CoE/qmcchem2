! Mu Jastrow
! --------------

! See Giner JCP 2021

BEGIN_PROVIDER [ double precision , jast_elec_Mu_value, (elec_num_8)  ]
implicit none
 BEGIN_DOC  
! J(i) = \sum_j a.rij/(1+b^2.rij) - \sum_A (a.riA/(1+a.riA))^2
! Eq (11) 
 END_DOC
 integer :: i,j
 double precision :: a, b, rij, tmp
 include '../constants.F'
 double precision :: mu
 mu = jast_mu_erf

 do i=1,elec_num
   jast_elec_Mu_value(i) =  jast_1b_value(i)
 enddo

 do j=1,elec_num 
  !DIR$ LOOP COUNT (50)
  do i=1,elec_num
    if(j==i)cycle
    rij = elec_dist(i,j)
    tmp = 0.5d0 * rij * (1.d0 - derf(mu*rij)) - 0.5d0/(dsqpi*mu) * dexp(-mu*mu*rij*rij)
    jast_elec_Mu_value(j) += 0.5d0*tmp ! symmetrization 
  enddo
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision , jast_elec_Mu_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Mu_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Mu_grad_z, (elec_num_8) ]
 implicit none
 BEGIN_DOC  
! Gradient of the Jastrow factor
! Eq (A1) 
 END_DOC

 integer :: i,j
 double precision :: a, b, rij, tmp, x, y, z
 include '../constants.F'
 double precision :: mu
 mu = jast_mu_erf

 do i=1,elec_num
  jast_elec_Mu_grad_x(i) = jast_1b_grad_x(i)
  jast_elec_Mu_grad_y(i) = jast_1b_grad_y(i)
  jast_elec_Mu_grad_z(i) = jast_1b_grad_z(i)
 enddo

 ! (grad of J(r12) with respect to xi, yi, zi)
 do i = 1, elec_num
  do j = 1, elec_num
   if(i==j)cycle
   rij = elec_dist(j,i)
   jast_elec_Mu_grad_x(i) += 0.5d0 * ( 1.d0 - derf(mu * rij) ) * elec_dist_inv(j,i) * (-1.d0) * elec_dist_vec_x(j,i)
   jast_elec_Mu_grad_y(i) += 0.5d0 * ( 1.d0 - derf(mu * rij) ) * elec_dist_inv(j,i) * (-1.d0) * elec_dist_vec_y(j,i)
   jast_elec_Mu_grad_z(i) += 0.5d0 * ( 1.d0 - derf(mu * rij) ) * elec_dist_inv(j,i) * (-1.d0) * elec_dist_vec_z(j,i)
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision , jast_elec_Mu_lapl, (elec_num_8) ]
 implicit none
 BEGIN_DOC  
! Laplacian of the Jastrow factor
! Eq (A10) 
 END_DOC

 integer :: i,j
 double precision :: a, b, rij, tmp, x, y, z
 include '../constants.F'
 double precision :: mu, x_ij, y_ij, z_ij, rij_inv
 mu = jast_mu_erf

 do i=1,elec_num
  jast_elec_Mu_lapl(i) = jast_1b_lapl(i)
 enddo

 do i=1, elec_num
  do j=1, elec_num
   if(j==i)cycle
   rij = elec_dist(j,i)
   rij_inv = elec_dist_inv(j,i)
   x_ij = elec_dist_vec_x(j,i)
   y_ij = elec_dist_vec_y(j,i)
   z_ij = elec_dist_vec_z(j,i)

   jast_elec_Mu_lapl(i) += (1.d0 - derf(mu*rij))*elec_dist_inv(j,i) - mu/dsqpi * dexp(-mu*mu*rij*rij)
  enddo
 enddo
END_PROVIDER


 BEGIN_PROVIDER [double precision, grad_j_mu_x,(elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_y,(elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_z,(elec_num, elec_num)]
 implicit none
 BEGIN_DOC
! Needed for 3-body terms
 END_DOC
 integer :: i,j
 double precision :: rij, mu,scal
 mu = jast_mu_erf
 grad_j_mu_x = 0.d0
 grad_j_mu_y = 0.d0
 grad_j_mu_z = 0.d0
 do j = 1, elec_num
  do i = 1, elec_num
   if(i==j)cycle
   rij = elec_dist(i,j)
   scal = 0.5d0 * ( 1.d0 - derf(mu * rij) ) * elec_dist_inv(i,j) 
   grad_j_mu_x(i,j) = (elec_coord_transp(1,i) - elec_coord_transp(1,j)) * scal
   grad_j_mu_y(i,j) = (elec_coord_transp(2,i) - elec_coord_transp(2,j)) * scal
   grad_j_mu_z(i,j) = (elec_coord_transp(3,i) - elec_coord_transp(3,j)) * scal
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, Energy_mu ]

  BEGIN_DOC
  ! E mu = < H_mu \Phi / \Phi >_{\Phi^2}
  END_DOC

  implicit none
  integer :: i

  double precision :: lapl
  lapl = 0.d0
  do i=1,elec_num
    lapl += psidet_grad_lapl(4,i)*psidet_inv + jast_elec_mu_lapl(i) + &
    2.d0*psidet_inv * (&
        psidet_grad_lapl(1,i)*jast_elec_mu_grad_x(i) +                   &
        psidet_grad_lapl(2,i)*jast_elec_mu_grad_y(i) +                   &
        psidet_grad_lapl(3,i)*jast_elec_mu_grad_z(i)  ) + ( &
        jast_elec_mu_grad_x(i)*jast_elec_mu_grad_x(i) + &
        jast_elec_mu_grad_y(i)*jast_elec_mu_grad_y(i) + &
        jast_elec_mu_grad_z(i)*jast_elec_mu_grad_z(i) )

  enddo
  Energy_mu = -0.5d0 * lapl + E_nucl + E_pot

END_PROVIDER

