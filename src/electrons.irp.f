use qmckl

BEGIN_PROVIDER [ double precision, xbrown, (elec_num_8,3) ]
  BEGIN_DOC
  ! Brownian step. Built in Brownian_step subroutine.
  END_DOC
  integer                        :: i,l
  xbrown = 0.d0
END_PROVIDER


 BEGIN_PROVIDER  [ integer, elec_alpha_num ]
&BEGIN_PROVIDER  [ integer, elec_alpha_num_8 ]
  implicit none
  BEGIN_DOC
  ! Number of alpha electrons
  END_DOC
  integer, external              :: mod_align
  elec_alpha_num = -1
  call get_electrons_elec_alpha_num(elec_alpha_num)
  if (elec_alpha_num <= 0) then
    call abrt(irp_here,'Number of alpha electrons should be > 0')
  endif
  elec_alpha_num_8 = mod_align(elec_alpha_num)

END_PROVIDER



 BEGIN_PROVIDER  [ integer, elec_beta_num ]
&BEGIN_PROVIDER  [ integer, elec_beta_num_8 ]
  implicit none
  BEGIN_DOC
  ! Number of beta electrons
  END_DOC
  integer, external              :: mod_align
  elec_beta_num = 0
  call get_electrons_elec_beta_num(elec_beta_num)
  if (elec_beta_num < 0) then
    call abrt(irp_here,'Number of beta electrons should be >= 0')
  endif
  elec_beta_num_8 = mod_align(elec_beta_num)

END_PROVIDER



 BEGIN_PROVIDER [ integer, elec_num ]
&BEGIN_PROVIDER [ integer, elec_num_8 ]
&BEGIN_PROVIDER [ integer, elec_num_1_8 ]
  implicit none
  BEGIN_DOC
  ! Number of electrons
  END_DOC
  integer, external              :: mod_align
  elec_num = elec_alpha_num + elec_beta_num
  ASSERT ( elec_num > 0 )
  elec_num_8 = mod_align(elec_num)
  elec_num_1_8 = mod_align(elec_num+1)

END_PROVIDER


BEGIN_PROVIDER   [ real, elec_coord_full, (elec_num+1,3,walk_num) ]
  implicit none
  BEGIN_DOC
  ! Electron coordinates of all walkers
  ! Component (elec_num+1,1,walk_num) contains the length realized by the walker.
  ! Initialized in init_walkers
  END_DOC
  integer :: ifirst_elec_coord_full
  common /common_elec_coord_full/ ifirst_elec_coord_full 
  integer                        :: i,k
  real, allocatable              :: buffer2(:,:,:)
  PROVIDE elec_coord_pool_size
  if ( is_worker ) then

    call get_elec_coord_full(elec_coord_full,size(elec_coord_full,1))

  else

    if (.not.do_prepare) then
      if (ifirst_elec_coord_full == 0) then
        allocate ( buffer2(elec_num+1,3,elec_coord_pool_size) )
        call get_electrons_elec_coord_pool(buffer2)
        do k=1,walk_num
          do i=1,elec_num+1
            elec_coord_full(i,1,k) = buffer2(i,1,k)
            elec_coord_full(i,2,k) = buffer2(i,2,k)
            elec_coord_full(i,3,k) = buffer2(i,3,k)
          enddo
        enddo
        deallocate ( buffer2 )
        ifirst_elec_coord_full=1
      endif
    else
      elec_coord_full = 0.
    endif

  endif
  double precision :: buffer(elec_num,3)
  integer(qmckl_exit_code) :: rc
  buffer(1:elec_num,1:3) = elec_coord_full(1:elec_num,1:3,1)
  rc = qmckl_set_electron_coord(qmckl_ctx, 'T', 1_8, buffer,  3_8*elec_num)
  call check_qmckl(rc, irp_here, qmckl_ctx)

END_PROVIDER

BEGIN_PROVIDER [ integer, elec_coord_pool_size ]
 implicit none
 BEGIN_DOC
 ! Size of the pool of electron coordinates
 END_DOC
 elec_coord_pool_size = walk_num_tot
 call get_electrons_elec_coord_pool_size(elec_coord_pool_size)
 call iinfo(irp_here,'elec_coord_pool_size',elec_coord_pool_size)

END_PROVIDER



BEGIN_PROVIDER  [ real, elec_coord, (elec_num_1_8,3) ]
  implicit none
  BEGIN_DOC
  ! Electron coordinates
  END_DOC
  integer                        :: i, k
  elec_coord = 0.
  do k=1,3
    do i=1,elec_num+1
      elec_coord(i,k) = elec_coord_full(i,k,walk_i)
    enddo
  enddo

  double precision :: buffer(elec_num,3)
  integer(qmckl_exit_code) :: rc
  buffer(1:elec_num,1:3) = elec_coord(1:elec_num,1:3)
  rc = qmckl_set_electron_coord(qmckl_ctx, 'T', 1_8, buffer,  3_8*elec_num)
  call check_qmckl(rc, irp_here, qmckl_ctx)

END_PROVIDER



BEGIN_PROVIDER [ real, elec_coord_transp, (8,elec_num)
  implicit none
  BEGIN_DOC
  ! Transposed array of elec_coord
  END_DOC
  integer                        :: i, k
  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    elec_coord_transp = 0.
  endif

  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i=1,elec_num
    elec_coord_transp(1,i) = elec_coord(i,1)
    elec_coord_transp(2,i) = elec_coord(i,2)
    elec_coord_transp(3,i) = elec_coord(i,3)
  enddo
END_PROVIDER



 BEGIN_PROVIDER [ real, elec_dist, (elec_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, elec_dist_vec_x, (elec_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, elec_dist_vec_y, (elec_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, elec_dist_vec_z, (elec_num_8,elec_num) ]
  implicit none
  BEGIN_DOC
  ! Electron-electron distances
  END_DOC
  integer                        :: ie1, ie2, l

  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    elec_dist = 0.
    !DIR$ VECTOR ALIGNED
    elec_dist_vec_x = 0.
    !DIR$ VECTOR ALIGNED
    elec_dist_vec_y = 0.
    !DIR$ VECTOR ALIGNED
    elec_dist_vec_z = 0.
  endif

  do ie2 = 1,elec_num
    real                           :: x, y, z
    real                           :: x2, y2, z2
    x = elec_coord(ie2,1)
    y = elec_coord(ie2,2)
    z = elec_coord(ie2,3)
    !DIR$ VECTOR ALIGNED
    do ie1 = 1,elec_num
      elec_dist_vec_x(ie1,ie2) = elec_coord(ie1,1) - x
      elec_dist_vec_y(ie1,ie2) = elec_coord(ie1,2) - y
      elec_dist_vec_z(ie1,ie2) = elec_coord(ie1,3) - z
      elec_dist(ie1,ie2) = sqrt(                                     &
          elec_dist_vec_x(ie1,ie2)*elec_dist_vec_x(ie1,ie2) +        &
          elec_dist_vec_y(ie1,ie2)*elec_dist_vec_y(ie1,ie2) +        &
          elec_dist_vec_z(ie1,ie2)*elec_dist_vec_z(ie1,ie2) )
    enddo
  enddo

END_PROVIDER



 BEGIN_PROVIDER [ real, nucl_elec_dist, (nucl_num_8,elec_num) ]
&BEGIN_PROVIDER [ real, nucl_elec_dist_vec, (3,nucl_num,elec_num) ]
  implicit none
  BEGIN_DOC
  !  Electron-nucleus distances  |r_elec - R_nucl|
  END_DOC
  integer                        :: i,j,l
  integer, save                  :: ifirst = 0

  if (ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    nucl_elec_dist = 0.
    !DIR$ VECTOR ALIGNED
    nucl_elec_dist_vec = 0.
  endif

  do i = 1,elec_num
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1,nucl_num
      nucl_elec_dist_vec(1,j,i) = elec_coord_transp(1,i) - nucl_coord(j,1)
      nucl_elec_dist_vec(2,j,i) = elec_coord_transp(2,i) - nucl_coord(j,2)
      nucl_elec_dist_vec(3,j,i) = elec_coord_transp(3,i) - nucl_coord(j,3)
    enddo
  enddo

  do i = 1,elec_num
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do j = 1,nucl_num
      nucl_elec_dist(j,i) = (elec_coord(i,1) - nucl_coord(j,1))      &
          * (elec_coord(i,1) - nucl_coord(j,1))                      &
          + (elec_coord(i,2) - nucl_coord(j,2))                      &
          * (elec_coord(i,2) - nucl_coord(j,2))                      &
          + (elec_coord(i,3) - nucl_coord(j,3))                      &
          * (elec_coord(i,3) - nucl_coord(j,3))
      nucl_elec_dist(j,i) = max(1.e-6,sqrt(nucl_elec_dist(j,i)))
    enddo
  enddo
END_PROVIDER



BEGIN_PROVIDER [ integer, elec_num_2, (2) ]

  BEGIN_DOC
  ! Number of alpha and beta electrons in an array
  END_DOC

  elec_num_2(1) = elec_alpha_num
  elec_num_2(2) = elec_beta_num

END_PROVIDER


BEGIN_PROVIDER [ integer, elec_spin, (elec_num) ]
  implicit none
  BEGIN_DOC
  ! Electron spin. +1 for alpha and -1 for beta
  END_DOC
  integer                        :: i
  do i=1,elec_alpha_num
    elec_spin(i) = 1
  enddo
  do i=elec_alpha_num+1,elec_num
    elec_spin(i) = -1
  enddo
END_PROVIDER



BEGIN_PROVIDER [ real, elec_dist_inv, (elec_num_8,elec_num) ]
  implicit none
  BEGIN_DOC
  ! 1/rij matrix
  END_DOC
  integer                        :: i,j
  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    elec_dist_inv = 0.
  endif
  do i=1,elec_num
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (200)
    do j=1,elec_num
      elec_dist_inv(j,i) = 1./(elec_dist(j,i)+1.e-12)
    enddo
    elec_dist_inv(i,i) = 0.
  enddo
END_PROVIDER


BEGIN_PROVIDER [ real, nucl_elec_dist_inv, (nucl_num_8,elec_num) ]
  implicit none
  BEGIN_DOC
  ! 1/rij matrix
  END_DOC
  integer                        :: i,j
  do j=1,elec_num
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do i=1,nucl_num
      nucl_elec_dist_inv(i,j) = 1./nucl_elec_dist(i,j)
    enddo
  enddo
END_PROVIDER


subroutine save_elec_coord_full
  implicit none
  BEGIN_DOC
! Save the electron coordinates to disk
  END_DOC
  integer                        :: i,k,l
  real, allocatable              :: buffer2(:,:,:)

  allocate ( buffer2(elec_num+1,3,elec_coord_pool_size) )
  k=0
  do l=1,elec_coord_pool_size
    k = k+1
    if (k == walk_num+1) then
      k=1
    endif
    do i=1,elec_num+1
      buffer2(i,1,l) = elec_coord_full(i,1,k)
      buffer2(i,2,l) = elec_coord_full(i,2,k)
      buffer2(i,3,l) = elec_coord_full(i,3,k)
    enddo
  enddo
  call ezfio_set_electrons_elec_coord_pool(buffer2)
  deallocate ( buffer2 )

end

