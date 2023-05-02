! Simple Jastrow
! --------------

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Simple_value, (elec_num_8)  ]
implicit none
 BEGIN_DOC  
! J(i) = \sum_j a.rij/(1+b.rij) - \sum_A (a.riA/(1+a.riA))^2
 END_DOC
 integer :: i,j
 double precision :: a, b, rij, tmp

 do i=1,elec_num
  jast_elec_Simple_value(i) = 0.d0
  !DIR$ LOOP COUNT (100)
  do j=1,nucl_num
    a = jast_pen(j)
    rij = nucl_elec_dist(j,i)
    tmp = a*rij/(1.d0+a*rij)
    jast_elec_Simple_value(i) -= tmp*tmp
  enddo
 enddo

 a = 0.5d0*jast_a_up_up
 b = jast_b_up_up

 ! parallele spin up-up
 do j=1,elec_alpha_num-1
  !DIR$ LOOP COUNT (50)
  do i=j+1,elec_alpha_num
    rij = elec_dist(i,j)
    tmp = a*rij/(1.d0+b*rij) 
    jast_elec_Simple_value(i) += tmp
    jast_elec_Simple_value(j) += tmp
  enddo
 enddo

 ! parallele spin dn-dn
 do j=elec_alpha_num+1,elec_num-1
  !DIR$ LOOP COUNT (50)
  do i=j+1,elec_num
    rij = elec_dist(i,j)
    tmp = a*rij/(1.d0+b*rij) 
    jast_elec_Simple_value(i) += tmp
    jast_elec_Simple_value(j) += tmp
  enddo
 enddo


 a = 0.5d0*jast_a_up_dn
 b = jast_b_up_dn

 ! anti-parallele spin
 do j=1,elec_alpha_num
  !DIR$ LOOP COUNT (50)
  do i=elec_alpha_num+1,elec_num
    rij = elec_dist(i,j)
    tmp = a*rij/(1.d0+b*rij) 
    jast_elec_Simple_value(i) += tmp
    jast_elec_Simple_value(j) += tmp
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision , jast_elec_Simple_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Simple_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Simple_grad_z, (elec_num_8) ]
 implicit none
 BEGIN_DOC  
! Gradient of the Jastrow factor
 END_DOC

 integer :: i,j
 double precision :: a, b, rij, tmp, x, y, z

 do i=1,elec_num
  jast_elec_Simple_grad_x(i) = 0.d0
  jast_elec_Simple_grad_y(i) = 0.d0
  jast_elec_Simple_grad_z(i) = 0.d0
  !DIR$ LOOP COUNT (100)
  do j=1,nucl_num
    a = jast_pen(j)
    rij = a*nucl_elec_dist(j,i)
    tmp = (a+a)*a/(1.d0+rij*(3.d0+rij*(3.d0+rij)))
    jast_elec_Simple_grad_x(i) -= nucl_elec_dist_vec(1,j,i)*tmp
    jast_elec_Simple_grad_y(i) -= nucl_elec_dist_vec(2,j,i)*tmp
    jast_elec_Simple_grad_z(i) -= nucl_elec_dist_vec(3,j,i)*tmp
  enddo
 enddo

 a = jast_a_up_up
 b = jast_b_up_up

 do j=1,elec_alpha_num
  !DIR$ LOOP COUNT (50)
  do i=1,j-1
    rij = elec_dist(i,j)
    tmp = a/(rij*(1.d0+b*rij*(2.d0+b*rij)))
    jast_elec_Simple_grad_x(i) += elec_dist_vec_x(i,j)*tmp
    jast_elec_Simple_grad_y(i) += elec_dist_vec_y(i,j)*tmp
    jast_elec_Simple_grad_z(i) += elec_dist_vec_z(i,j)*tmp
    jast_elec_Simple_grad_x(j) -= elec_dist_vec_x(i,j)*tmp
    jast_elec_Simple_grad_y(j) -= elec_dist_vec_y(i,j)*tmp
    jast_elec_Simple_grad_z(j) -= elec_dist_vec_z(i,j)*tmp
  enddo
 enddo

 do j=elec_alpha_num+1,elec_num
  !DIR$ LOOP COUNT (50)
  do i=elec_alpha_num+1,j-1
    rij = elec_dist(i,j)
    tmp = a/(rij*(1.d0+b*rij*(2.d0+b*rij)))
    jast_elec_Simple_grad_x(i) += elec_dist_vec_x(i,j)*tmp
    jast_elec_Simple_grad_y(i) += elec_dist_vec_y(i,j)*tmp
    jast_elec_Simple_grad_z(i) += elec_dist_vec_z(i,j)*tmp
    jast_elec_Simple_grad_x(j) -= elec_dist_vec_x(i,j)*tmp
    jast_elec_Simple_grad_y(j) -= elec_dist_vec_y(i,j)*tmp
    jast_elec_Simple_grad_z(j) -= elec_dist_vec_z(i,j)*tmp
  enddo
 enddo

 a = jast_a_up_dn
 b = jast_b_up_dn

 do j=1,elec_alpha_num
  !DIR$ LOOP COUNT (50)
  do i=elec_alpha_num+1,elec_num
    rij = elec_dist(i,j)
    tmp = a/(rij*(1.d0+b*rij*(2.d0+b*rij)))
    jast_elec_Simple_grad_x(i) += elec_dist_vec_x(i,j)*tmp
    jast_elec_Simple_grad_y(i) += elec_dist_vec_y(i,j)*tmp
    jast_elec_Simple_grad_z(i) += elec_dist_vec_z(i,j)*tmp
    jast_elec_Simple_grad_x(j) -= elec_dist_vec_x(i,j)*tmp
    jast_elec_Simple_grad_y(j) -= elec_dist_vec_y(i,j)*tmp
    jast_elec_Simple_grad_z(j) -= elec_dist_vec_z(i,j)*tmp
  enddo
 enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Simple_lapl, (elec_num_8) ]

  BEGIN_DOC  
  ! Laplacian of the Jastrow factor
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: a, b, rij, tmp

  do i = 1, elec_num
    jast_elec_Simple_lapl(i) = 0.d0
    !DIR$ LOOP COUNT (100)
    do j=1,nucl_num
      a = jast_pen(j)
      rij = a*nucl_elec_dist(j,i)
      tmp = 6.d0*a*a/(1.d0+rij*(4.d0+rij*(6.d0+rij*(4.d0+rij))))
      jast_elec_Simple_lapl(i) -= tmp
    enddo
  enddo

  a = jast_a_up_up + jast_a_up_up
  b = jast_b_up_up
  do j = 1, elec_alpha_num
    !DIR$ LOOP COUNT (50)
    do i = 1, j-1
      rij = b*elec_dist(i,j)
      tmp = a/(elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
      jast_elec_Simple_lapl(i) += tmp
      jast_elec_Simple_lapl(j) += tmp
    enddo
  enddo
  do j = elec_alpha_num+1, elec_num
  !DIR$ LOOP COUNT (100)
  do i = j+1, elec_num
    rij = b*elec_dist(i,j)
    tmp = a/(elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
    jast_elec_Simple_lapl(i) += tmp
    jast_elec_Simple_lapl(j) += tmp
  enddo
 enddo

  a = jast_a_up_dn + jast_a_up_dn
  b = jast_b_up_dn
  do j = 1, elec_alpha_num
    !DIR$ LOOP COUNT (100)
    do i = elec_alpha_num+1, elec_num
      rij = b*elec_dist(i,j)
      tmp = a/(elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
      jast_elec_Simple_lapl(i) += tmp
      jast_elec_Simple_lapl(j) += tmp
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER[ double precision, jast_elec_Simple_deriv_nucPar, (nucl_num) ]
  implicit none
  BEGIN_DOC  
  ! Variation of the Jastrow factor with respect to nuclear parameters
  END_DOC

  integer          :: i, j
  double precision :: a, rij, tmp1, tmp2

  do j = 1, nucl_num
    a    = jast_pen(j)
    tmp2 = 0.d0
    !DIR$ LOOP COUNT (100)
    do i = 1, elec_num
      rij  = nucl_elec_dist(j,i)
      tmp1 = (1.d0+a*rij)*(1.d0+a*rij)*(1.d0+a*rij)
      tmp2 += rij*rij/tmp1
    end do
    jast_elec_Simple_deriv_nucPar(j) = -2.d0 * a * tmp2
  end do

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision , Jsimple_value ]
&BEGIN_PROVIDER [ double precision , Jsimple_value_inv ]

  implicit none
  integer          :: i
  double precision :: argexpo

  integer, save    :: ifirst = 0
  double precision :: dshift = 0.d0

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Simple_value(i)
  enddo

  if (ifirst == 0) then
    dshift = argexpo
    ifirst = 1
  endif
  argexpo -= dshift

  Jsimple_value     = dexp(argexpo)
  Jsimple_value_inv = 1.d0 / Jsimple_value

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_elec_Simple_lapl_reg, (elec_num_8) ]

  BEGIN_DOC  
  ! Regularized laplacian of the Simple Jastrow:
  !      1 / [ r (1 + br)^3 ] - 1 / r = -[ 3 + 3br + (br)^2 ] b / (1 + br)^3 
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: a, b, rij, tmp

  ! 1-body term
  do i = 1, elec_num
    jast_elec_Simple_lapl_reg(i) = 0.d0
    !DIR$ LOOP COUNT (100)
    do j = 1, nucl_num
      a = jast_pen(j)
      rij = a * nucl_elec_dist(j,i)
      tmp = 6.d0*a*a / (1.d0+rij*(4.d0+rij*(6.d0+rij*(4.d0+rij))))
      jast_elec_Simple_lapl_reg(i) -= tmp
    enddo
  enddo

  ! 2-body term
  b = jast_b_up_up
  do i = 1, elec_num
    !DIR$ LOOP COUNT (50)
    do j = 1, elec_num
      if(i==j) cycle
      rij = b * elec_dist(j,i)
      a   = 3.d0 + rij * ( 3.d0 + rij ) 
      tmp = a * b / ( 1.d0 + rij * a )
      jast_elec_Simple_lapl_reg(i) -= tmp
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jsimple_nonh ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= ( psidet_right_grad_lapl(1,i) * jast_elec_Simple_grad_x(i) &
           + psidet_right_grad_lapl(2,i) * jast_elec_Simple_grad_y(i) &
           + psidet_right_grad_lapl(3,i) * jast_elec_Simple_grad_z(i) ) * psidet_right_inv
  enddo
  deltaE_Jsimple_nonh = tmp 

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jsimple_grad ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Simple_grad_x(i) * jast_elec_Simple_grad_x(i) &
         + jast_elec_Simple_grad_y(i) * jast_elec_Simple_grad_y(i) &
         + jast_elec_Simple_grad_z(i) * jast_elec_Simple_grad_z(i)
  enddo
  deltaE_Jsimple_grad = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jsimple_lapl ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Simple_lapl(i)
  enddo
  deltaE_Jsimple_lapl = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jsimple_lapl_reg ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Simple_lapl_reg(i)
  enddo
  deltaE_Jsimple_lapl_reg = 0.5d0 * tmp

END_PROVIDER

! ---
 
 
