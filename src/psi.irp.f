
! --- 

BEGIN_PROVIDER [ double precision, psi_value ]

  BEGIN_DOC
  ! Value of the wave function
  END_DOC
  
  implicit none

  if(use_lr) then

    psi_value = coef_psi_right * psidet_right_value * jast_value &
              + coef_psi_left  * psidet_left_value  * jast_value_inv

  else

    psi_value = psidet_right_value * jast_value

  endif

  if(psi_value == 0.d0) then
    call abrt(irp_here,"Value of the wave function is 0.")
  endif

END_PROVIDER

! --- 

BEGIN_PROVIDER [ double precision, psi_value_inv ]
  implicit none
  BEGIN_DOC
! 1./psi_value
  END_DOC
  psi_value_inv = 1.d0/psi_value
END_PROVIDER

! --- 

BEGIN_PROVIDER [ double precision, psi_value_inv2 ]
 implicit none
 BEGIN_DOC
 ! 1./(psi_value)**2
 END_DOC
 psi_value_inv2 = psi_value_inv*psi_value_inv
END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, psi_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, psi_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, psi_grad_z, (elec_num_8) ]

  BEGIN_DOC
  !
  ! grad[a Phi_R e^{J} + b Phi_L e^{-J}] = a [grad[Phi_R] + Phi_R grad[J]] e^{+J}  
  !                                      + b [grad[Phi_L] - Phi_L grad[J]] e^{-J}
  !
  ! grad[Phi e^{sgn_jast x J}] = [grad[Phi] + Phi grad[sgn_jast x J]] e^{sgn_jast x J}
  !
  END_DOC
  
  implicit none
  integer          :: j
  double precision :: gradr_x, gradr_y, gradr_z, gradl_x, gradl_y, gradl_z

  if(use_lr) then

    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num

      gradr_x = (psidet_right_grad_lapl(1,j) + psidet_right_value * jast_grad_jast_inv_x(j)) * jast_value
      gradr_y = (psidet_right_grad_lapl(2,j) + psidet_right_value * jast_grad_jast_inv_y(j)) * jast_value
      gradr_z = (psidet_right_grad_lapl(3,j) + psidet_right_value * jast_grad_jast_inv_z(j)) * jast_value

      gradl_x = (psidet_left_grad_lapl(1,j) - psidet_left_value * jast_grad_jast_inv_x(j)) * jast_value_inv
      gradl_y = (psidet_left_grad_lapl(2,j) - psidet_left_value * jast_grad_jast_inv_y(j)) * jast_value_inv
      gradl_z = (psidet_left_grad_lapl(3,j) - psidet_left_value * jast_grad_jast_inv_z(j)) * jast_value_inv

      psi_grad_x(j) = coef_psi_right * gradr_x + coef_psi_left * gradl_x
      psi_grad_y(j) = coef_psi_right * gradr_y + coef_psi_left * gradl_y
      psi_grad_z(j) = coef_psi_right * gradr_z + coef_psi_left * gradl_z

    enddo

  else
    ! here, sgn_jast is in jast_grad_jast_inv_x, jast_grad_jast_inv_y, jast_grad_jast_inv_z
    ! and jast_value = e{sgn_jast x J}
    
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num
      psi_grad_x(j) = (psidet_right_grad_lapl(1,j) + psidet_right_value * jast_grad_jast_inv_x(j)) * jast_value
      psi_grad_y(j) = (psidet_right_grad_lapl(2,j) + psidet_right_value * jast_grad_jast_inv_y(j)) * jast_value
      psi_grad_z(j) = (psidet_right_grad_lapl(3,j) + psidet_right_value * jast_grad_jast_inv_z(j)) * jast_value
    enddo

  endif
  
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, psi_lapl, (elec_num_8) ]

  BEGIN_DOC
  ! Laplacian of the wave function
  END_DOC
  
  implicit none
  integer          :: i, j
  double precision :: jlapl1, jlapl2
  double precision :: jgrad_x, jgrad_y, jgrad_z, tmp_r, tmp_l

  if(use_lr) then

    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num

      jlapl1  = jast_lapl1(j)
      jlapl2  = jast_lapl2(j)
      jgrad_x = jast_grad_jast_inv_x(j)
      jgrad_y = jast_grad_jast_inv_y(j)
      jgrad_z = jast_grad_jast_inv_z(j)

      tmp_r = psidet_right_grad_lapl(4,j)                     &
            + psidet_right_value * (jlapl1 + jlapl2)          &
            + 2.d0 * ( psidet_right_grad_lapl(1,j) * jgrad_x  &
                     + psidet_right_grad_lapl(2,j) * jgrad_y  &
                     + psidet_right_grad_lapl(3,j) * jgrad_z )

      tmp_l = psidet_left_grad_lapl(4,j)                     &
            + psidet_left_value * (-jlapl1 + jlapl2)         &
            - 2.d0 * ( psidet_left_grad_lapl(1,j) * jgrad_x  &
                     + psidet_left_grad_lapl(2,j) * jgrad_y  &
                     + psidet_left_grad_lapl(3,j) * jgrad_z )

      psi_lapl(j) = coef_psi_right * tmp_r * jast_value + coef_psi_left * tmp_l * jast_value_inv

    enddo

  else

    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num
      psi_lapl(j) = jast_value * ( psidet_right_grad_lapl(4,j)                                    &
                                 + psidet_right_value * jast_lapl_jast_inv(j)                     &
                                 + 2.d0 * ( psidet_right_grad_lapl(1,j) * jast_grad_jast_inv_x(j) &
                                          + psidet_right_grad_lapl(2,j) * jast_grad_jast_inv_y(j) &
                                          + psidet_right_grad_lapl(3,j) * jast_grad_jast_inv_z(j) ))
    enddo

  endif
  
END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, psi_grad_psi_inv_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, psi_grad_psi_inv_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, psi_grad_psi_inv_z, (elec_num_8) ]

  BEGIN_DOC
  ! grad(psi)/psi
  END_DOC
  
  implicit none
  integer :: j

  if(use_lr) then

    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num
      psi_grad_psi_inv_x(j) = psi_grad_x(j) * psi_value_inv 
      psi_grad_psi_inv_y(j) = psi_grad_y(j) * psi_value_inv 
      psi_grad_psi_inv_z(j) = psi_grad_z(j) * psi_value_inv 
    enddo

  else
    
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1, elec_num
      psi_grad_psi_inv_x(j) = psidet_right_grad_lapl(1,j) * psidet_right_inv + jast_grad_jast_inv_x(j)
      psi_grad_psi_inv_y(j) = psidet_right_grad_lapl(2,j) * psidet_right_inv + jast_grad_jast_inv_y(j)
      psi_grad_psi_inv_z(j) = psidet_right_grad_lapl(3,j) * psidet_right_inv + jast_grad_jast_inv_z(j)
    enddo

  endif
  
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, psi_lapl_psi_inv, (elec_num_8) ]

  BEGIN_DOC
  ! (Laplacian psi) / psi
  END_DOC
  
  implicit none
  integer :: i, j

  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (100)
  do j = 1, elec_num
    psi_lapl_psi_inv(j) = psi_lapl(j) * psi_value_inv
  enddo

END_PROVIDER

! ---

