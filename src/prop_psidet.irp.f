
! ---

BEGIN_PROVIDER [ double precision, E_kin_elec_psidet_right, (elec_num) ]

  BEGIN_DOC
  ! Electronic Kinetic energy of the right determinantal part only: -1/2 (Lapl.Psidet)/Psidet
  END_DOC

  implicit none
  integer :: i

  do i = 1, elec_num
    E_kin_elec_psidet_right(i) = -0.5d0 * psidet_right_grad_lapl(4,i) * psidet_right_inv
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, E_kin_psidet_right ]

  BEGIN_DOC
  ! Electronic Kinetic energy of the right determinantal part only: -1/2 (Lapl.Psidet)/Psidet
  END_DOC

  implicit none
  integer :: i

  E_kin_psidet_right = 0.d0
  do i = 1, elec_num
    E_kin_psidet_right -= 0.5d0 * psidet_right_grad_lapl(4,i) * psidet_right_inv
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, E_kin_elec_psidet_left, (elec_num) ]

  BEGIN_DOC
  ! Electronic Kinetic energy of the left determinantal part only: -1/2 (Lapl.Psidet)/Psidet
  END_DOC

  implicit none
  integer :: i

  do i = 1, elec_num
    E_kin_elec_psidet_left(i) = -0.5d0 * psidet_left_grad_lapl(4,i) * psidet_left_inv
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, E_kin_psidet_left ]

  BEGIN_DOC
  ! Electronic Kinetic energy of the left determinantal part only: -1/2 (Lapl.Psidet)/Psidet
  END_DOC

  implicit none
  integer :: i

  E_kin_psidet_left = 0.d0
  do i = 1, elec_num
    E_kin_psidet_left -= 0.5d0 * psidet_left_grad_lapl(4,i) * psidet_left_inv
  enddo

END_PROVIDER

! ---

