! Mu Jastrow + 1body
! -------------------

! ---

 BEGIN_PROVIDER [ double precision, jast_Mu_1b_value ]
&BEGIN_PROVIDER [ double precision, jast_Mu_1b_value_inv ]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Mu_1b_value(i)
  enddo

  jast_Mu_1b_value     = dexp(sgn_jast*argexpo)
  jast_Mu_1b_value_inv = 1.d0 / jast_Mu_1b_value

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Mu_1b_value, (elec_num_8)  ]

  implicit none
  integer :: i

  do i = 1, elec_num
    jast_elec_Mu_1b_value(i) = jast_elec_Mu_value(i) + jast_1b_value(i)
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision , jast_elec_Mu_1b_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Mu_1b_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Mu_1b_grad_z, (elec_num_8) ]

  implicit none
  integer :: i

  do i = 1, elec_num
    jast_elec_Mu_1b_grad_x(i) = jast_elec_Mu_grad_x(i) + jast_1b_grad_x(i)
    jast_elec_Mu_1b_grad_y(i) = jast_elec_Mu_grad_y(i) + jast_1b_grad_y(i)
    jast_elec_Mu_1b_grad_z(i) = jast_elec_Mu_grad_z(i) + jast_1b_grad_z(i)
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Mu_1b_lapl, (elec_num_8) ]

  implicit none
  integer :: i

  do i = 1, elec_num
    jast_elec_Mu_1b_lapl(i) = jast_elec_Mu_lapl(i) + jast_1b_lapl(i)
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Mu_1b_lapl_reg, (elec_num_8) ]

  implicit none
  integer :: i

  do i = 1, elec_num
    jast_elec_Mu_1b_lapl_reg(i) = jast_elec_Mu_lapl_reg(i) + jast_1b_lapl(i)
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_lapl ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_1b_lapl(i)
  enddo
  deltaE_Jmu1b_lapl = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_lapl_reg ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_1b_lapl_reg(i)
  enddo
  deltaE_Jmu1b_lapl_reg = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_nonh ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= ( psidet_right_grad_lapl(1,i) * jast_elec_Mu_1b_grad_x(i) &
           + psidet_right_grad_lapl(2,i) * jast_elec_Mu_1b_grad_y(i) &
           + psidet_right_grad_lapl(3,i) * jast_elec_Mu_1b_grad_z(i) ) * psidet_right_inv
  enddo
  deltaE_Jmu1b_nonh = tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_mJmu1b_nonh ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= ( psidet_left_grad_lapl(1,i) * jast_elec_Mu_1b_grad_x(i) &
           + psidet_left_grad_lapl(2,i) * jast_elec_Mu_1b_grad_y(i) &
           + psidet_left_grad_lapl(3,i) * jast_elec_Mu_1b_grad_z(i) ) * psidet_left_inv
  enddo
  deltaE_mJmu1b_nonh = tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_grad ]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_1b_grad_x(i) * jast_elec_Mu_1b_grad_x(i) &
         + jast_elec_Mu_1b_grad_y(i) * jast_elec_Mu_1b_grad_y(i) &
         + jast_elec_Mu_1b_grad_z(i) * jast_elec_Mu_1b_grad_z(i)
  enddo
  deltaE_Jmu1b_grad = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_grad_no3b ]

  implicit none
  integer          :: i, j
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    !DIR$ LOOP COUNT(100)
    do j = 1, elec_num
      if(i==j) cycle
      ! + sign for i <--> j
      tmp += grad_j_mu_x(j,i) * ( jast_1b_grad_x(i) - jast_1b_grad_x(j) ) &
           + grad_j_mu_y(j,i) * ( jast_1b_grad_y(i) - jast_1b_grad_y(j) ) &
           + grad_j_mu_z(j,i) * ( jast_1b_grad_z(i) - jast_1b_grad_z(j) )
    enddo
  enddo

  ! x 0.5 for double counting
  deltaE_Jmu1b_grad_no3b = deltaE_Jmu_grad_2b + deltaE_J1b_grad + 0.5d0 * tmp 

END_PROVIDER

! ---

