 BEGIN_PROVIDER [ integer, ao_num ]
&BEGIN_PROVIDER [ integer, ao_num_8 ]
  implicit none
  BEGIN_DOC
! Number of atomic orbitals
  END_DOC
  integer, external              :: mod_align

  ao_num = -1
  call get_ao_basis_ao_num(ao_num)
  if (ao_num <= 0) then
    call abrt(irp_here,'Number of contracted gaussians should be > 0')
  endif
  call iinfo(irp_here,'ao_num',ao_num)
  ao_num_8 = mod_align(ao_num)
END_PROVIDER


BEGIN_PROVIDER [ integer, ao_prim_num, (ao_num_8) ]
  implicit none
  BEGIN_DOC
! Number of primitives per atomic orbital
  END_DOC

  ao_prim_num = 0
  call get_ao_basis_ao_prim_num(ao_prim_num)
  integer                        :: i
  character*(80)                 :: message
  do i=1,ao_num
    if (ao_prim_num(i) <= 0) then
      write(message,'(A,I6,A)') 'Number of primitives of contraction ',i,' should be > 0'
      call abrt(irp_here,message)
    endif
  enddo
  call iset_order(ao_prim_num,ao_nucl_sort_idx,ao_num)
  call iinfo(irp_here,'sum_ao_prim_num',sum(ao_prim_num))

END_PROVIDER


 BEGIN_PROVIDER [ integer, ao_nucl, (ao_num_8) ]
&BEGIN_PROVIDER [ integer, ao_nucl_sort_idx, (ao_num_8) ]
&BEGIN_PROVIDER [ integer, ao_nucl_idx, (2,nucl_num) ]
  implicit none
  BEGIN_DOC
! Nucleus on which the atomic orbital is centered
  END_DOC

  ao_nucl = -1
  call get_ao_basis_ao_nucl(ao_nucl)

  character*(80)                 :: message
  character*(30)                 :: range
  integer                        :: i,k

  do i=1,ao_num
    if ( (ao_nucl(i) <= 0) .or. (ao_nucl(i) > nucl_num) ) then
      write(range,'(A,I5,A)') '(1,',nucl_num,')'
      write(message,'(A,I6,A)') 'Contraction ',i,' should be centered on a nucleus in the range'//trim(range)
      call abrt(irp_here,message)
    endif
  enddo
  do i=1,ao_num_8
    ao_nucl_sort_idx(i) = i
  enddo

  call insertion_isort(ao_nucl,ao_nucl_sort_idx,ao_num)
  ao_nucl_idx(1,ao_nucl(1)) = 1
  ao_nucl_idx(2,ao_nucl(1)) = 1
  do i=2,ao_num
    k = ao_nucl(i)
    if (k == ao_nucl(i-1)) then
      ao_nucl_idx(2,k) = i
    else
      ao_nucl_idx(1,k) = i
      ao_nucl_idx(2,k) = i
    endif
  enddo

END_PROVIDER


 BEGIN_PROVIDER [ integer, ao_power, (ao_num,4) ]
&BEGIN_PROVIDER [ logical, ao_power_is_zero, (ao_num) ]
  implicit none
  BEGIN_DOC
! x,y,z powers of the atomic orbital
!
! 4 contains the sum of powers
!
! ao_power_is_zero is true where ao_power(:,4) == 0
  END_DOC
  ao_power = 0
  call get_ao_basis_ao_power(ao_power)

  character*(80)                 :: message
  integer                        :: i,j
  do i=1,3
    do j=1,ao_num
      if (ao_power(j,i) < 0) then
        write(message,'(A,I1,A,I6,A)') 'Power ',i,' of contraction ',j,' should be > 0'
        call abrt(irp_here,message)
      endif
    enddo
    call iset_order(ao_power(1,i),ao_nucl_sort_idx,ao_num)
  enddo
  do j=1,ao_num
    ao_power(j,4) = ao_power(j,1)+ao_power(j,2)+ao_power(j,3)
    ao_power_is_zero(j) = ao_power(j,4) == 0
  enddo

END_PROVIDER


BEGIN_PROVIDER [ integer, ao_power_transp, (4,ao_num) ]
  implicit none
  BEGIN_DOC
! Transposed ao_power
  END_DOC
  integer                        :: i,j
  do i=1,ao_num
    do j=1,4
      ao_power_transp(j,i) = ao_power(i,j)
    enddo
  enddo
END_PROVIDER



BEGIN_PROVIDER [ integer , ao_power_max ]
  implicit none
  BEGIN_DOC
! Maximum power among x, y and z
  END_DOC
  ao_power_max = maxval(ao_power_max_nucl)
END_PROVIDER


BEGIN_PROVIDER [ integer , ao_power_max_nucl, (nucl_num_8,4) ]
  implicit none
  BEGIN_DOC
! Maximum powers of x, y and z per nucleus
  END_DOC
  integer                        :: i, j
  ao_power_max_nucl = 0

  integer                        :: inucl
  do j=1,3
    do i=1,ao_num
      inucl = ao_nucl(i)
      ao_power_max_nucl(inucl,j) = max(ao_power(i,j),ao_power_max_nucl(inucl,j))
    enddo
  enddo
  do inucl=1,nucl_num
    ao_power_max_nucl(inucl,4) = max(ao_power_max_nucl(inucl,1),     &
        max(ao_power_max_nucl(inucl,2),ao_power_max_nucl(inucl,3)) )
  enddo
END_PROVIDER


 BEGIN_PROVIDER [ integer, ao_prim_num_max ]
&BEGIN_PROVIDER [ integer, ao_prim_num_max_2 ]
&BEGIN_PROVIDER [ integer, ao_prim_num_max_8 ]
&BEGIN_PROVIDER [ integer, ao_prim_num_max_pow2 ]
  implicit none
  BEGIN_DOC
! Max Number of primitives per atomic orbital
  END_DOC
  integer, external              :: mod_align

  ao_prim_num_max = maxval(ao_prim_num)
  call iinfo(irp_here,'ao_prim_num_max',ao_prim_num_max)
  ao_prim_num_max_8 = mod_align(ao_prim_num_max)
  ao_prim_num_max_2 = mod(ao_prim_num_max,2)
  if (ao_prim_num_max_2 == 0) then
    ao_prim_num_max_2 = ao_prim_num_max
  else
    ao_prim_num_max_2 = ao_prim_num_max + 2 - ao_prim_num_max_2
  endif
  ao_prim_num_max_pow2 = 16
  do while (ao_prim_num_max_pow2 < ao_prim_num_max)
    ao_prim_num_max_pow2 *= 2
  enddo

END_PROVIDER


 BEGIN_PROVIDER [ real, ao_expo, (ao_num_8,ao_prim_num_max) ]
&BEGIN_PROVIDER [ real, ao_coef, (ao_num_8,ao_prim_num_max) ]
  implicit none
  BEGIN_DOC
! Exponents and coefficients of the atomic orbitals.
! AO coefficients are such that the AOs are normalized, and the primitives
! are also normalized
  END_DOC

  ao_expo = 0.
  real, allocatable              :: buf(:,:)
  allocate (buf(ao_num,ao_prim_num_max))
  buf = 0.
  call get_ao_basis_ao_expo(buf)
  do j=1,ao_prim_num_max
    do i=1,ao_num
      ao_expo(i,j) = buf(i,j)
    enddo
    call set_order(ao_expo(1,j),ao_nucl_sort_idx,ao_num)
  enddo

  integer                        :: i,j
  do i=1,ao_num
    do j=1,ao_prim_num(i)
      if (ao_expo(i,j) <= 0.) then
        character*(80)                 :: message
        write(message,'(A,I6,A,I6,A)') 'Exponent ',j,' of contracted gaussian ',i,' is < 0'
        call abrt(irp_here,message)
      endif
    enddo
  enddo

  ao_coef = 0.
  buf = 0.
  call get_ao_basis_ao_coef(buf)
  do j=1,ao_prim_num_max
    do i=1,ao_num
      ao_coef(i,j) = buf(i,j)
    enddo
    call set_order(ao_coef(1,j),ao_nucl_sort_idx,ao_num)
  enddo
  deallocate(buf)

  real,allocatable               :: buffer(:)
  integer, allocatable           :: order(:)
  allocate (buffer(ao_prim_num_max),order(ao_prim_num_max))
  do i=1,ao_num
    do j=1,ao_prim_num(i)
      buffer(j) = ao_expo(i,j)
      order(j) = j
    enddo
    call sort(buffer,order,ao_prim_num(i))
    do j=1,ao_prim_num(i)
      ao_expo(i,j) = buffer(j)
      buffer(j) = ao_coef(i,j)
    enddo
    do j=1,ao_prim_num(i)
      ao_coef(i,order(j)) = buffer(j)
    enddo
  enddo
  deallocate(buffer,order)

! Normalization of the AO coefficients
! ------------------------------------
  real                           :: norm, norm2
  real                           :: goverlap
  integer                        :: pow(3), l
  do i=1,ao_num
    do j=1,ao_prim_num(i)
      pow(1) = ao_power_transp(1,i)
      pow(2) = ao_power_transp(2,i)
      pow(3) = ao_power_transp(3,i)
      norm = goverlap(ao_expo(i,j),ao_expo(i,j),pow)
      ao_coef(i,j) = ao_coef(i,j)/sqrt(norm)
    enddo
  enddo

END_PROVIDER


 BEGIN_PROVIDER [ real, ao_expo_unique, (2*sum(ao_prim_num)) ]
&BEGIN_PROVIDER [ integer, ao_expo_unique_nucl_idx, (2,nucl_num) ]
&BEGIN_PROVIDER [ integer, ao_expo_unique_idx_1, (ao_num_8) ]
&BEGIN_PROVIDER [ integer, ao_expo_unique_idx, (ao_prim_num_max_pow2,ao_num) ]
  implicit none
  BEGIN_DOC
! Exponents and coefficients of the atomic orbitals
  END_DOC

  integer                        :: k, kk, kkstart, kk2, i, j
  integer, allocatable           :: order(:)
  allocate ( order(2*sum(ao_prim_num)) )

  ao_expo_unique_idx = 1
  kk=1
  do k=1,nucl_num
    kkstart = kk
    do i=1,ao_num
      if (ao_nucl(i) == k) then
        do j=1,ao_prim_num(i)
          order(kk) = kk-kkstart+1
          ao_expo_unique(kk) = ao_expo(i,j)
          kk += 1
        enddo
      endif
    enddo
    call sort(ao_expo_unique(kkstart),order(kkstart),kk-kkstart)
    order(kk) = 0
    kk = kkstart
    kk2 = kk+1
    do while (order(kk2) /= 0)
      if ( ao_expo_unique(kk2) /= ao_expo_unique(kk) ) then
        kk += 1
        ao_expo_unique(kk) = ao_expo_unique(kk2)
      endif
      kk2+=1
    enddo
    ao_expo_unique_nucl_idx(1,k) = kkstart
    ao_expo_unique_nucl_idx(2,k) = kk
    kk += 1
  enddo
  deallocate(order)

  ao_expo_unique_idx = 0
  do i=1,ao_num
    do j=1,ao_prim_num(i)
      do k=ao_expo_unique_nucl_idx(1,ao_nucl(i)),ao_expo_unique_nucl_idx(2,ao_nucl(i))
        if (ao_expo_unique(k) == ao_expo(i,j)) then
          ao_expo_unique_idx(j,i) = k
          exit
        endif
      enddo
    enddo
    ao_expo_unique_idx_1(i) = ao_expo_unique_idx(1,i)
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_coef_transp, (ao_prim_num_max_8,ao_num) ]
  implicit none
  BEGIN_DOC
! Transposed of the ao_coef matrix
  END_DOC
  integer                        :: i,j
  integer, save                  :: ifirst = 0
  if (ifirst==0) then
    ifirst = 1
    ao_coef_transp = 0.
  endif
  call transpose(ao_coef,ao_num_8,ao_coef_transp,ao_prim_num_max_8,ao_num,ao_prim_num_max)
END_PROVIDER


BEGIN_PROVIDER [ real, ao_expo_transp, (ao_prim_num_max_8,ao_num) ]
  implicit none
  BEGIN_DOC
! Transposed of the ao_expo matrix
  END_DOC
  integer                        :: i,j
  integer, save                  :: ifirst = 0
  if (ifirst==0) then
    ifirst = 1
    ao_expo_transp = 0.
  endif
  call transpose(ao_expo,ao_num_8,ao_expo_transp,ao_prim_num_max_8,ao_num,ao_prim_num_max)
END_PROVIDER

