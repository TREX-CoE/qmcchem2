
! ---

 BEGIN_PROVIDER [ double precision, jpsi_value ]
&BEGIN_PROVIDER [ double precision, jpsi_value_inv ]

  BEGIN_DOC  
  ! Value of the Jastrow factor
  END_DOC

  include '../types.F'

  implicit none
  integer          :: i
  integer, save    :: ifirst = 0
  double precision :: dshift = 0.d0
  double precision :: argexpo

  jpsi_value = 1.d0

  if(do_jpsi) then

BEGIN_TEMPLATE
    if(jpsi_type == t_$X) then
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

    jpsi_value = dexp(argexpo)

  endif

 ASSERT (jpsi_value > 0.d0)
 jpsi_value_inv = 1.d0 / jpsi_value

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, jpsi_grad_jpsi_inv_x, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jpsi_grad_jpsi_inv_y, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jpsi_grad_jpsi_inv_z, (elec_num) ]
&BEGIN_PROVIDER [ double precision, jpsi_grad_x,          (elec_num) ]
&BEGIN_PROVIDER [ double precision, jpsi_grad_y,          (elec_num) ]
&BEGIN_PROVIDER [ double precision, jpsi_grad_z,          (elec_num) ]

  BEGIN_DOC  
  ! Grad(J)/J 
  END_DOC

  include '../types.F'

  implicit none
  integer       :: i, l
  integer, save :: ifirst = 0

  if(ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (200)
    do i = 1, elec_num
      jpsi_grad_jpsi_inv_x(i) = 0.d0
      jpsi_grad_jpsi_inv_y(i) = 0.d0
      jpsi_grad_jpsi_inv_z(i) = 0.d0
      jpsi_grad_x(i)          = 0.d0
      jpsi_grad_y(i)          = 0.d0
      jpsi_grad_z(i)          = 0.d0
    enddo
  endif

  if(do_jpsi) then

BEGIN_TEMPLATE
    if(jpsi_type == t_$X) then
      !DIR$ VECTOR ALIGNED
      !DIR$ LOOP COUNT (200)
      do i = 1, elec_num
        jpsi_grad_jpsi_inv_x(i) = jast_elec_$X_grad_x(i)
        jpsi_grad_jpsi_inv_y(i) = jast_elec_$X_grad_y(i)
        jpsi_grad_jpsi_inv_z(i) = jast_elec_$X_grad_z(i)
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
      jpsi_grad_x(i) = jpsi_grad_jpsi_inv_x(i) * jpsi_value
      jpsi_grad_y(i) = jpsi_grad_jpsi_inv_y(i) * jpsi_value
      jpsi_grad_z(i) = jpsi_grad_jpsi_inv_z(i) * jpsi_value
    enddo

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, jpsi_lapl,          (elec_num) ]
&BEGIN_PROVIDER [ double precision, jpsi_lapl_jpsi_inv, (elec_num) ]

  BEGIN_DOC  
  ! Lapl(J)/J
  END_DOC

  include '../types.F'
 
  implicit none
  integer       :: i 
  integer, save :: ifirst = 0

  if(ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    jpsi_lapl_jpsi_inv = 0.d0
    jpsi_lapl          = 0.d0
  endif
 
  if(do_jpsi) then

BEGIN_TEMPLATE
    if(jpsi_type == t_$X) then
      !DIR$ VECTOR ALIGNED
      !DIR$ LOOP COUNT (200)
      do i = 1, elec_num
        jpsi_lapl_jpsi_inv(i) = jast_elec_$X_lapl(i)                              &
                              + jpsi_grad_jpsi_inv_x(i) * jpsi_grad_jpsi_inv_x(i) &
                              + jpsi_grad_jpsi_inv_y(i) * jpsi_grad_jpsi_inv_y(i) &
                              + jpsi_grad_jpsi_inv_z(i) * jpsi_grad_jpsi_inv_z(i)
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

    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (256)
    do i = 1, elec_num
      jpsi_lapl(i) = jpsi_lapl_jpsi_inv(i) * jpsi_value
    enddo

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_elec_Jpsi_lapl, (elec_num) ]

  BEGIN_DOC  
  ! Lapl(J)
  END_DOC

  include '../types.F'
 
  implicit none
  integer       :: i 
  integer, save :: ifirst = 0

  if(ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    jast_elec_Jpsi_lapl = 0.d0
  endif
 
  if(do_jpsi) then

BEGIN_TEMPLATE
    if(jpsi_type == t_$X) then
      !DIR$ VECTOR ALIGNED
      !DIR$ LOOP COUNT (200)
      do i = 1, elec_num
        jast_elec_Jpsi_lapl(i) = jast_elec_$X_lapl(i)
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

!   print *, jast_elec_Jpsi_lapl(1) , jast_elec_Qmckl_lapl(1)
  endif

END_PROVIDER


! ---

BEGIN_PROVIDER [ double precision, deltaE_Jpsi_lapl ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Jpsi_lapl(i)
  enddo
  deltaE_Jpsi_lapl = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jpsi_nonh ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= ( psidet_right_grad_lapl(1,i) * jpsi_grad_jpsi_inv_x(i) &
           + psidet_right_grad_lapl(2,i) * jpsi_grad_jpsi_inv_y(i) &
           + psidet_right_grad_lapl(3,i) * jpsi_grad_jpsi_inv_z(i) ) * psidet_right_inv
  enddo
  deltaE_Jpsi_nonh = tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_mJpsi_nonh ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= ( psidet_left_grad_lapl(1,i) * jpsi_grad_jpsi_inv_x(i) &
           + psidet_left_grad_lapl(2,i) * jpsi_grad_jpsi_inv_y(i) &
           + psidet_left_grad_lapl(3,i) * jpsi_grad_jpsi_inv_z(i) ) * psidet_left_inv
  enddo
  deltaE_mJpsi_nonh = tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jpsi_grad ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jpsi_grad_jpsi_inv_x(i) * jpsi_grad_jpsi_inv_x(i) &
         + jpsi_grad_jpsi_inv_y(i) * jpsi_grad_jpsi_inv_y(i) &
         + jpsi_grad_jpsi_inv_z(i) * jpsi_grad_jpsi_inv_z(i)
  enddo
  deltaE_Jpsi_grad = 0.5d0 * tmp

END_PROVIDER

! ---


