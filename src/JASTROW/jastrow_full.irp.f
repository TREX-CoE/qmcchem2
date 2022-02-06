 BEGIN_PROVIDER [ double precision, jast_value ]
&BEGIN_PROVIDER [ double precision, jast_value_inv ]
 implicit none
 include '../types.F'
 BEGIN_DOC  
! Value of the Jastrow factor
 END_DOC

 integer, save :: ifirst = 0
 integer :: i
 double precision :: dshift = 0.d0

 jast_value = 1.d0
 if (do_jast) then
   double precision :: argexpo
BEGIN_TEMPLATE
   if (jast_type == t_$X) then
     argexpo = 0.d0
     !DIR$ VECTOR ALIGNED
     !DIR$ LOOP COUNT (200)
     do i=1,elec_num
       argexpo += jast_elec_$X_value(i)
     enddo
 !   argexpo = argexpo/dble(elec_num)
   endif
SUBST [X]
Simple ;;
Core   ;;
Mu     ;;
END_TEMPLATE
   if (ifirst == 0) then
     dshift = argexpo
     ifirst = 1
   endif
   argexpo -= dshift
   jast_value = exp(argexpo)
 endif
 ASSERT (jast_value > 0.d0)
 jast_value_inv = 1.d0/jast_value

END_PROVIDER

 BEGIN_PROVIDER [ double precision, jast_grad_jast_inv_x, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_grad_jast_inv_y, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_grad_jast_inv_z, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_grad_x, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_grad_y, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_grad_z, (elec_num) ]
 implicit none
 include '../types.F'
 BEGIN_DOC  
! Grad(J)/J 
 END_DOC

 integer :: i,l
 integer, save :: ifirst = 0

 if (ifirst == 0) then
   ifirst = 1
   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (200)
   do i=1,elec_num
    jast_grad_jast_inv_x(i) = 0.d0
    jast_grad_jast_inv_y(i) = 0.d0
    jast_grad_jast_inv_z(i) = 0.d0
    jast_grad_x(i) = 0.d0
    jast_grad_y(i) = 0.d0
    jast_grad_z(i) = 0.d0
   enddo
 endif

 if (do_jast) then

BEGIN_TEMPLATE
   if ( jast_type == t_$X ) then
      !DIR$ VECTOR ALIGNED
      !DIR$ LOOP COUNT (200)
      do i=1,elec_num
       jast_grad_jast_inv_x(i) = jast_elec_$X_grad_x(i)
       jast_grad_jast_inv_y(i) = jast_elec_$X_grad_y(i)
       jast_grad_jast_inv_z(i) = jast_elec_$X_grad_z(i)
      enddo
   endif
SUBST [ X ]
Simple ;;
Core   ;;
Mu     ;;
END_TEMPLATE
   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (200)
   do i=1,elec_num
    jast_grad_x(i) = jast_grad_jast_inv_x(i)*jast_value
    jast_grad_y(i) = jast_grad_jast_inv_y(i)*jast_value
    jast_grad_z(i) = jast_grad_jast_inv_z(i)*jast_value
   enddo

 endif
END_PROVIDER


 BEGIN_PROVIDER [ double precision, jast_lapl, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jast_lapl_jast_inv, (elec_num) ]
 implicit none
 include '../types.F'
 BEGIN_DOC  
! Lapl(J)/J
 END_DOC

 integer :: i
 integer,save :: ifirst = 0

 if (ifirst == 0) then
   ifirst = 1
   !DIR$ VECTOR ALIGNED
   jast_lapl_jast_inv = 0.d0
 endif

 if (do_jast) then

BEGIN_TEMPLATE
   if (jast_type == t_$X) then
     !DIR$ VECTOR ALIGNED
     !DIR$ LOOP COUNT (200)
     do i=1,elec_num
       jast_lapl_jast_inv(i) = (jast_elec_$X_lapl(i) + &
            jast_grad_jast_inv_x(i)*jast_grad_jast_inv_x(i) + &
            jast_grad_jast_inv_y(i)*jast_grad_jast_inv_y(i) + &
            jast_grad_jast_inv_z(i)*jast_grad_jast_inv_z(i) )
     enddo
   endif
SUBST [X]
Simple ;;
Core ;;
Mu   ;;
END_TEMPLATE

   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (256)
   do i=1,elec_num
     jast_lapl(i) = jast_lapl_jast_inv(i) * jast_value
   enddo

 endif
END_PROVIDER

