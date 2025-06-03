! Core Jastrow
! --------------

 BEGIN_PROVIDER [ double precision, jast_elec_Core_expo, (nucl_num) ]
&BEGIN_PROVIDER [ double precision, jast_elec_Core_range, (nucl_num) ]
 implicit none
 BEGIN_DOC
! Exponent of the core jastrow factor per nucleus
 END_DOC
 integer :: i
 do i=1,nucl_num
   if (nucl_charge(i) < 2.5d0) then
     jast_elec_Core_expo(i) = 0.d0
     jast_elec_Core_range(i) = 0.d0
   else
     double precision :: rc
     double precision, parameter :: thresh = 0.5d0  ! function = thresh at rc
     rc = min(0.8d0,max(4.0d0/nucl_charge(i), 0.25d0))
     jast_elec_Core_expo(i) = -1.d0/rc**2 * log(thresh)
     jast_elec_Core_range(i) = dsqrt(15.d0/jast_elec_Core_expo(i))
   endif
   call dinfo(irp_here, 'expo', jast_elec_Core_expo(i))
 enddo

END_PROVIDER




BEGIN_PROVIDER [ double precision , jast_elec_Core_value, (elec_num)  ]
implicit none
 BEGIN_DOC
! J(i) = \sum_j a.rij/(1+b^2.rij) - \sum_A (a.riA/(1+a.riA))^2
 END_DOC
 integer :: i,j,k
 double precision :: a, b, rij, tmp, a2
 double precision :: f1

 jast_elec_Core_value = 0.d0

 do k=1,nucl_num
  if (jast_elec_Core_range(k) == 0.d0) then
    cycle
  endif

  a = 0.5d0
  a2 = jast_core_a1(k)
  b  = jast_core_b1(k)

  do j=1,elec_alpha_num
   if (nucl_elec_dist(k,j) > jast_elec_Core_range(k)) then
     cycle
   endif
   do i=elec_alpha_num+1,elec_num
     if (nucl_elec_dist(k,i) > jast_elec_Core_range(k)) then
      cycle
     endif
     rij = elec_dist(i,j)
     f1 = exp(-jast_elec_Core_expo(k)*(nucl_elec_dist(k,i)*nucl_elec_dist(k,i)+nucl_elec_dist(k,j)*nucl_elec_dist(k,j)))
     tmp = f1*(a*rij/(1.d0+b*rij) - a2)
     jast_elec_Core_value(i) = jast_elec_Core_value(i) + 0.5d0*tmp
     jast_elec_Core_value(j) = jast_elec_Core_value(j) + 0.5d0*tmp
   enddo
  enddo

 enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision , jast_elec_Core_grad_x, (elec_num) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Core_grad_y, (elec_num) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Core_grad_z, (elec_num) ]
 implicit none
 BEGIN_DOC
! Gradient of the Jastrow factor
 END_DOC

 integer :: i,j,k
 double precision :: a, b, rij, tmp, x, y, z, f1, a2

 jast_elec_Core_grad_x = 0.d0
 jast_elec_Core_grad_y = 0.d0
 jast_elec_Core_grad_z = 0.d0


 do k=1,nucl_num
   if (jast_elec_Core_range(k) == 0.d0) then
    cycle
   endif

   a = 0.5d0
   a2 = jast_core_a1(k)
   b  = jast_core_b1(k)

   do j=1,elec_alpha_num
    if (nucl_elec_dist(k,j) > jast_elec_Core_range(k)) then
     cycle
    endif
    do i=elec_alpha_num+1,elec_num
      if (nucl_elec_dist(k,i) > jast_elec_Core_range(k)) then
       cycle
      endif
      rij = elec_dist(i,j)
      f1 = exp(-jast_elec_Core_expo(k)*(nucl_elec_dist(k,i)*nucl_elec_dist(k,i)+nucl_elec_dist(k,j)*nucl_elec_dist(k,j)))
      tmp = f1*(a/(rij*(1.d0+b*rij*(2.d0+b*rij))))
      f1 = -2.d0*jast_elec_Core_expo(k)* f1* (a*rij/(1.d0+b*rij) -a2)
      jast_elec_Core_grad_x(i) +=  elec_dist_vec_x(i,j)*tmp + nucl_elec_dist_vec(1,k,i)*f1
      jast_elec_Core_grad_y(i) +=  elec_dist_vec_y(i,j)*tmp + nucl_elec_dist_vec(2,k,i)*f1
      jast_elec_Core_grad_z(i) +=  elec_dist_vec_z(i,j)*tmp + nucl_elec_dist_vec(3,k,i)*f1
      jast_elec_Core_grad_x(j) += -elec_dist_vec_x(i,j)*tmp + nucl_elec_dist_vec(1,k,j)*f1
      jast_elec_Core_grad_y(j) += -elec_dist_vec_y(i,j)*tmp + nucl_elec_dist_vec(2,k,j)*f1
      jast_elec_Core_grad_z(j) += -elec_dist_vec_z(i,j)*tmp + nucl_elec_dist_vec(3,k,j)*f1
    enddo
   enddo

 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision , jast_elec_Core_lapl, (elec_num) ]
 implicit none
 BEGIN_DOC
! Laplacian of the Jastrow factor
 END_DOC

 integer :: i,j,k
 double precision :: a, b, rij, tmp, x, y, z,f1, alpha, a2

 jast_elec_Core_lapl = 0.d0

 do k=1,nucl_num
   if (jast_elec_Core_range(k) == 0.d0) then
    cycle
   endif

   a = 0.5d0
   a2 = jast_core_a1(k)
   b  = jast_core_b1(k)

   do j=1,elec_alpha_num
    if (nucl_elec_dist(k,j) > jast_elec_Core_range(k)) then
     cycle
    endif
    do i=elec_alpha_num+1,elec_num
      if (nucl_elec_dist(k,i) > jast_elec_Core_range(k)) then
       cycle
      endif
      f1 = exp(-jast_elec_Core_expo(k)*(nucl_elec_dist(k,i)*nucl_elec_dist(k,i)+nucl_elec_dist(k,j)*nucl_elec_dist(k,j)))
      rij = b*elec_dist(i,j)
      tmp = (a+a)/(elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
      jast_elec_Core_lapl(i) += tmp*f1
      jast_elec_Core_lapl(j) += tmp*f1

      rij = elec_dist(i,j)
      tmp = f1* ( a*rij/(1.d0+b*rij) -a2 )
      jast_elec_Core_lapl(i) += tmp*( 4.d0*(nucl_elec_dist(k,i)*jast_elec_Core_expo(k))**2 &
            -6.d0*jast_elec_Core_expo(k))
      jast_elec_Core_lapl(j) += tmp*( 4.d0*(nucl_elec_dist(k,j)*jast_elec_Core_expo(k))**2 &
            -6.d0*jast_elec_Core_expo(k))


      tmp =  4.d0*jast_elec_Core_expo(k)*f1*(a/(rij*(1.d0+b*rij*(2.d0+b*rij))) )
      jast_elec_Core_lapl(i) -= tmp*(nucl_elec_dist_vec(1,k,i)*elec_dist_vec_x(i,j) &
                                   + nucl_elec_dist_vec(2,k,i)*elec_dist_vec_y(i,j) &
                                   + nucl_elec_dist_vec(3,k,i)*elec_dist_vec_z(i,j))
      jast_elec_Core_lapl(j) += tmp*(nucl_elec_dist_vec(1,k,j)*elec_dist_vec_x(i,j) &
                                   + nucl_elec_dist_vec(2,k,j)*elec_dist_vec_y(i,j) &
                                   + nucl_elec_dist_vec(3,k,j)*elec_dist_vec_z(i,j))
    enddo
   enddo

 enddo

END_PROVIDER
