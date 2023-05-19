
subroutine rho_hf_val_der_lap(r, rho_val, rho_der, rho_lap)

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: rho_val, rho_der(3), rho_lap
  integer                       :: i
  real                          :: r_tmp(3)
  real                          :: mo_val_tmp, mo_der_tmp(3), mo_lap_tmp
  double precision              :: mo_val, mo_der(3), mo_lap, mo_der_square
  double precision              :: dx, dy, dz, r2

  !dx         = r(1)
  !dy         = r(2)
  !dz         = r(3)
  !r2         = dx*dx + dy*dy + dz*dz
  !rho_val    = exp(-r2)
  !rho_der(1) = -2.d0 * rho_val * dx
  !rho_der(2) = -2.d0 * rho_val * dy
  !rho_der(3) = -2.d0 * rho_val * dz
  !rho_lap    = (4.d0 * r2 - 6.d0) * rho_val
  !return

  rho_val    = 0.d0
  rho_der(1) = 0.d0
  rho_der(2) = 0.d0
  rho_der(3) = 0.d0
  rho_lap    = 0.d0

  do i = 1, elec_beta_num
    !call get_mo_val_der_lap(i, r, mo_val, mo_der, mo_lap)
    r_tmp(1) = real(r(1))
    r_tmp(2) = real(r(2))
    r_tmp(3) = real(r(3))
    call get_mo_val_der_lap(i, r_tmp, mo_val_tmp, mo_der_tmp, mo_lap_tmp)
    mo_val    = dble(mo_val_tmp   )
    mo_der(1) = dble(mo_der_tmp(1))
    mo_der(2) = dble(mo_der_tmp(2))
    mo_der(3) = dble(mo_der_tmp(3))
    mo_lap    = dble(mo_lap_tmp   )

    rho_val       = rho_val + mo_val * mo_val
    rho_der(1)    = rho_der(1) + 2.d0 * mo_val * mo_der(1)
    rho_der(2)    = rho_der(2) + 2.d0 * mo_val * mo_der(2)
    rho_der(3)    = rho_der(3) + 2.d0 * mo_val * mo_der(3)
    mo_der_square = mo_der(1) * mo_der(1) + mo_der(2) * mo_der(2) + mo_der(3) * mo_der(3)
    rho_lap       = rho_lap + 2.d0 * (mo_der_square + mo_val * mo_lap)
  enddo

  rho_val    = 2.d0 * rho_val   
  rho_der(1) = 2.d0 * rho_der(1)
  rho_der(2) = 2.d0 * rho_der(2)
  rho_der(3) = 2.d0 * rho_der(3)
  rho_lap    = 2.d0 * rho_lap   

  if(elec_alpha_num .gt. elec_beta_num) then
    do i = elec_beta_num+1, elec_alpha_num
      !call get_mo_val_der_lap(i, r, mo_val, mo_der, mo_lap)
      r_tmp(1) = real(r(1))
      r_tmp(2) = real(r(2))
      r_tmp(3) = real(r(3))
      call get_mo_val_der_lap(i, r_tmp, mo_val_tmp, mo_der_tmp, mo_lap_tmp)
      mo_val    = dble(mo_val_tmp   )
      mo_der(1) = dble(mo_der_tmp(1))
      mo_der(2) = dble(mo_der_tmp(2))
      mo_der(3) = dble(mo_der_tmp(3))
      mo_lap    = dble(mo_lap_tmp   )

      rho_val       = rho_val + mo_val * mo_val
      rho_der(1)    = rho_der(1) + 2.d0 * mo_val * mo_der(1)
      rho_der(2)    = rho_der(2) + 2.d0 * mo_val * mo_der(2)
      rho_der(3)    = rho_der(3) + 2.d0 * mo_val * mo_der(3)
      mo_der_square = mo_der(1) * mo_der(1) + mo_der(2) * mo_der(2) + mo_der(3) * mo_der(3)
      rho_lap       = rho_lap + 2.d0 * (mo_der_square + mo_val * mo_lap)
    enddo
  endif

  !print*, rho_val
  !print*, rho_der
  !print*, rho_lap

  return
end subroutine rho_hf_val_der_lap

! ---

subroutine get_mo_val_der_lap(i, r, mo_val, mo_der, mo_lap)

  BEGIN_DOC
  !
  ! i is the MO (1-->mo_num)
  !
  END_DOC

  implicit none
  integer, intent(in)  :: i
  real,    intent(in)  :: r(3)
  real,    intent(out) :: mo_val, mo_der(3), mo_lap
  integer              :: mu
  real                 :: ao_val, ao_der(3), ao_lap

  PROVIDE mo_coef_aux

  mo_val    = 0. 
  mo_der(1) = 0. 
  mo_der(2) = 0. 
  mo_der(3) = 0. 
  mo_lap    = 0. 

  do mu = 1, ao_num

    call get_ao_val_der_lap(mu, r, ao_val, ao_der, ao_lap)

    mo_val    = mo_val    + mo_coef_aux(mu,i) * ao_val
    mo_der(1) = mo_der(1) + mo_coef_aux(mu,i) * ao_der(1)
    mo_der(2) = mo_der(2) + mo_coef_aux(mu,i) * ao_der(2)
    mo_der(3) = mo_der(3) + mo_coef_aux(mu,i) * ao_der(3)
    mo_lap    = mo_lap    + mo_coef_aux(mu,i) * ao_lap

  enddo

  return

end subroutine get_mo_val_der_lap

! ---

BEGIN_PROVIDER [real, mo_coef_aux_input, (ao_num_8,mo_tot_num)]

  implicit none
  integer           :: i, j
  real, allocatable :: buffer(:,:)

  allocate(buffer(ao_num,mo_tot_num))
  buffer = 0.

  call get_mo_basis_mo_coef_aux(buffer)

  do i = 1, mo_tot_num
    do j = 1, ao_num
      mo_coef_aux_input(j,i) = buffer(j,i)
    enddo
    call set_order(mo_coef_aux_input(1,i), ao_nucl_sort_idx, ao_num)
    do j = ao_num+1, ao_num_8
      mo_coef_aux_input(j,i) = 0.
    enddo
  enddo

  deallocate(buffer)

END_PROVIDER

! ---

BEGIN_PROVIDER [real, mo_coef_aux, (ao_num_8,mo_num_8)]

  implicit none
  integer :: i, j

  do j = 1, mo_num
    do i = 1, ao_num_8
      mo_coef_aux(i,j) = mo_coef_aux_input(i,j)
    enddo
  enddo
  do j = mo_num+1, mo_num_8
    !DIR$ VECTOR ALIGNED
    do i = 1, ao_num_8
      mo_coef(i,j) = 0.
    enddo
  enddo

  ! Input MOs are not needed any more
  FREE mo_coef_aux_input

END_PROVIDER

! ---

