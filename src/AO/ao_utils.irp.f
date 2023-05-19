
subroutine get_ao_val_der_lap(i, r, ao_val, ao_der, ao_lap)

  BEGIN_DOC
  !
  ! i is the AO (1-->ao_num)
  !
  END_DOC

  implicit none
  integer, intent(in)  :: i
  real,    intent(in)  :: r(3)
  real,    intent(out) :: ao_val, ao_der(3), ao_lap
  integer              :: l, inuc
  real                 :: d(3)
  real                 :: ao_axis_val, ao_axis_der(3), ao_axis_lap
  real                 :: ao_1d_val  , ao_1d_der(3)  , ao_1d_lap

  inuc = ao_nucl(i)

  !print *, 'i        = ', i
  !print *, 'inuc     = ', inuc
  !print *, 'position = ', r

  d(1) = r(1) - nucl_coord(inuc,1)
  d(2) = r(2) - nucl_coord(inuc,2)
  d(3) = r(3) - nucl_coord(inuc,3)

  call get_ao_useful(i, d, ao_axis_val, ao_axis_der, ao_axis_lap, ao_1d_val, ao_1d_der, ao_1d_lap)

  ao_val = ao_1d_val * ao_axis_val

  do l = 1, 3
    ao_der(l) = ao_1d_val * ao_axis_der(l) + ao_1d_der(l) * ao_axis_val
  enddo

  ao_lap = ao_1d_val * ao_axis_lap + ao_1d_lap * ao_axis_val
  do l = 1, 3
    ao_lap = ao_lap + 2. * ao_1d_der(l) * ao_axis_der(l)
  enddo

  !print*, ao_val
  !print*, ao_der
  !print*, ao_lap

  return
end subroutine get_ao_val_der_lap

! ---

subroutine get_ao_useful(i, d, ao_axis_val, ao_axis_der, ao_axis_lap, ao_1d_val, ao_1d_der, ao_1d_lap)

  include '../types.F'

  implicit none
  integer, intent(in)  :: i
  real,    intent(in)  :: d(3)
  real,    intent(out) :: ao_axis_val, ao_axis_der(3), ao_axis_lap
  real,    intent(out) :: ao_1d_val  , ao_1d_der(3)  , ao_1d_lap
  integer              :: k, l, ii, jj, kk, inuc
  real                 :: r2, ao_1d_prim_expo(ao_prim_num_max), ao_1d_prim_der(ao_prim_num_max,3), ao_1d_prim_lap(ao_prim_num_max)
  real                 :: ao_axis_power(-2:ao_power_max,3)
  real                 :: real_of_int1(-1:10), real_of_int2(-2:10)

  inuc = ao_nucl(i)
  r2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)

  ! ---

  ! useful exponential terms
  ao_1d_prim_expo = 0.
  do k = 1, ao_prim_num(i)
    ao_1d_prim_expo(k) = exp(-ao_expo(i,k)*r2)
  enddo

  ao_1d_val = 0.
  do k = 1, ao_prim_num_max
    ao_1d_val = ao_1d_val + ao_coef(i,k) * ao_1d_prim_expo(k)
  enddo

  do l = 1, 3
    do k = 1, ao_prim_num_max
      ao_1d_prim_der(k,l) = -2. * d(l) * ao_expo(i,k) * ao_1d_prim_expo(k)
    enddo
  enddo
  do l = 1, 3
    ao_1d_der(l) = 0.
    do k = 1, ao_prim_num_max
      ao_1d_der(l) = ao_1d_der(l) + ao_coef(i,k) * ao_1d_prim_der(k,l)
    enddo
  enddo

  do k = 1, ao_prim_num_max
    ao_1d_prim_lap(k) = ao_1d_prim_expo(k) * ao_expo(i,k) * (4. * ao_expo(i,k) * r2 - 6.)
  enddo
  ao_1d_lap = 0.
  do k = 1, ao_prim_num_max
    ao_1d_lap = ao_1d_lap + ao_coef(i,k) * ao_1d_prim_lap(k)
  enddo

  ! ---

  ! useful powers
  do l = 1, 3
    ao_axis_power(-2,l) = 0.
    ao_axis_power(-1,l) = 0.
    ao_axis_power( 0,l) = 0.
    ao_axis_power( 0,l) = 1.
    do k = 1, ao_power_max_nucl(inuc,l)
      ao_axis_power(k,l) = d(l) * ao_axis_power(k-1,l)
    enddo
  enddo

  ii = ao_power_transp(1,i)
  jj = ao_power_transp(2,i)
  kk = ao_power_transp(3,i)

  ! Cartesian polynomial part of the atomic orbitals.
  ao_axis_val = ao_axis_power(ii,1) * ao_axis_power(jj,2) * ao_axis_power(kk,3)

  data real_of_int1 /0.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
  ! Gradients of the cartesian polynomial part of the atomic orbitals.
  ao_axis_der(1) = real_of_int1(ii) * ao_axis_power(ii-1,1) * ao_axis_power(jj  ,2) * ao_axis_power(kk  ,3)
  ao_axis_der(2) = real_of_int1(jj) * ao_axis_power(ii  ,1) * ao_axis_power(jj-1,2) * ao_axis_power(kk  ,3)
  ao_axis_der(3) = real_of_int1(kk) * ao_axis_power(ii  ,1) * ao_axis_power(jj  ,2) * ao_axis_power(kk-1,3)

  data real_of_int2 /0.,0.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
  ! Laplacian of the cartesian polynomial part of the atomic orbitals.
  ao_axis_lap = real_of_int2(ii) * real_of_int2(ii-1) * ao_axis_power(ii-2,1) * ao_axis_power(jj  ,2) * ao_axis_power(kk  ,3) &
              + real_of_int2(jj) * real_of_int2(jj-1) * ao_axis_power(ii  ,1) * ao_axis_power(jj-2,2) * ao_axis_power(kk  ,3) &
              + real_of_int2(kk) * real_of_int2(kk-1) * ao_axis_power(ii  ,1) * ao_axis_power(jj  ,2) * ao_axis_power(kk-2,3)

  ! ---

  return
end subroutine get_ao_useful

