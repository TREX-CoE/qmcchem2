BEGIN_PROVIDER [ real, mo_value_p, (mo_tot_num) ]
  implicit none

  BEGIN_DOC
  ! Values of the molecular orbitals
  END_DOC

  integer                        :: i, j, k

  do j=1,mo_num
    mo_value_p(j) = 0.d0
  enddo
  do k=1,ao_num
    do j=1,mo_num
      mo_value_p(j) = mo_value_p(j)+mo_coef_transp(j,k)*ao_value_p(k)
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ real, mo_left_value_p, (mo_tot_num) ]

  BEGIN_DOC
  ! Values of the left molecular orbitals
  END_DOC

  implicit none
  integer :: i, j, k

  do j = 1, mo_num
    mo_left_value_p(j) = 0.d0
  enddo
  do k = 1, ao_num
    do j = 1, mo_num
      mo_left_value_p(j) = mo_left_value_p(j) + mo_left_coef_transp(j,k) * ao_value_p(k)
    enddo
  enddo

END_PROVIDER

! ---

