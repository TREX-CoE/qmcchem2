 BEGIN_PROVIDER [ double precision, jast_value ]
&BEGIN_PROVIDER [ double precision, jast_value_inv ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Value of the Jastrow factor
 END_DOC

 integer          :: i
 integer, save    :: ifirst = 0
 double precision :: dshift = 0.d0
 double precision :: argexpo

 jast_value = 1.d0
 if (do_jast) then

BEGIN_TEMPLATE
     if (jast_type == t_$X) then
       argexpo = 0.d0
       !DIR$ VECTOR ALIGNED
       !DIR$ LOOP COUNT (200)
       do i = 1, elec_num
         argexpo += jast_elec_$X_value(i)
       enddo
     endif
SUBST [X]
Core   ;;
Mu     ;;
Mu_1b  ;;
Muenv  ;;
Mur    ;;
Qmckl  ;;
Simple ;;
END_TEMPLATE

     if (ifirst == 0) then
       dshift = argexpo
       ifirst = 1
     endif

     jast_value = dexp(sgn_jast*argexpo)
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
        do i = 1, elec_num
          jast_grad_jast_inv_x(i) = sgn_jast * jast_elec_$X_grad_x(i)
          jast_grad_jast_inv_y(i) = sgn_jast * jast_elec_$X_grad_y(i)
          jast_grad_jast_inv_z(i) = sgn_jast * jast_elec_$X_grad_z(i)
        enddo
     endif
SUBST [ X ]
Core   ;;
Mu     ;;
Mu_1b  ;;
Muenv  ;;
Mur    ;;
Qmckl  ;;
Simple ;;
END_TEMPLATE

   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (200)
   do i = 1, elec_num
     jast_grad_x(i) = jast_grad_jast_inv_x(i) * jast_value
     jast_grad_y(i) = jast_grad_jast_inv_y(i) * jast_value
     jast_grad_z(i) = jast_grad_jast_inv_z(i) * jast_value
   enddo

 endif
END_PROVIDER


! ---

BEGIN_PROVIDER [ double precision, jast_lapl_jast_inv, (elec_num) ]

  BEGIN_DOC
  !
  ! Lapl[e^{sgn_jast x J}] / e^{sgn_jast x J} = sgn_jast * Lapl[J] + grad(sgn_jast x J) . grad(sgn_jast x J)
  !
  END_DOC

  include '../types.F'

  implicit none

  integer       :: i
  integer, save :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    jast_lapl_jast_inv = 0.d0
  endif

  if(do_jast) then

BEGIN_TEMPLATE
      if(jast_type == t_$X) then
        !DIR$ VECTOR ALIGNED
        !DIR$ LOOP COUNT (200)
        do i = 1, elec_num
          jast_lapl_jast_inv(i) = sgn_jast * jast_elec_$X_lapl(i)                   &
                                + jast_grad_jast_inv_x(i) * jast_grad_jast_inv_x(i) &
                                + jast_grad_jast_inv_y(i) * jast_grad_jast_inv_y(i) &
                                + jast_grad_jast_inv_z(i) * jast_grad_jast_inv_z(i)
        enddo
      endif
SUBST [X]
Core   ;;
Mu     ;;
Mu_1b  ;;
Muenv  ;;
Mur    ;;
Qmckl  ;;
Simple ;;
END_TEMPLATE

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_lapl1, (elec_num) ]

  BEGIN_DOC
  !
  ! Lapl[J]
  !
  END_DOC

  include '../types.F'

  implicit none

  integer       :: i
  integer, save :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    jast_lapl1 = 0.d0
  endif

  if(do_jast) then

BEGIN_TEMPLATE
    if(jast_type == t_$X) then
      !DIR$ VECTOR ALIGNED
      !DIR$ LOOP COUNT (200)
      do i = 1, elec_num
        jast_lapl1(i) = jast_elec_$X_lapl(i)
      enddo
    endif
SUBST [X]
Core   ;;
Mu     ;;
Mu_1b  ;;
Muenv  ;;
Mur    ;;
Qmckl  ;;
Simple ;;
END_TEMPLATE

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_lapl2, (elec_num) ]

  BEGIN_DOC
  !
  ! grad(J) . grad(J)
  !
  END_DOC

  include '../types.F'

  implicit none

  integer       :: i
  integer, save :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    jast_lapl2 = 0.d0
  endif

  if(do_jast) then

    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (200)
    do i = 1, elec_num
      jast_lapl2(i) = jast_grad_jast_inv_x(i) * jast_grad_jast_inv_x(i) &
                    + jast_grad_jast_inv_y(i) * jast_grad_jast_inv_y(i) &
                    + jast_grad_jast_inv_z(i) * jast_grad_jast_inv_z(i)
    enddo

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_lapl, (elec_num) ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Lapl(J)
 END_DOC

 integer :: i

 jast_lapl(:) = jast_lapl_jast_inv(:) * jast_value

END_PROVIDER

