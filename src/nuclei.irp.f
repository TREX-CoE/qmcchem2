 BEGIN_PROVIDER [ integer, nucl_num ]
&BEGIN_PROVIDER [ integer, nucl_num_8 ]
  implicit none
  BEGIN_DOC
! Number of nuclei
  END_DOC

  nucl_num = -1
  call get_nuclei_nucl_num(nucl_num)
  if (nucl_num <= 0) then
    call abrt(irp_here,'Number of nuclei should be > 0')
  endif
  integer, external              :: mod_align
  nucl_num_8 = mod_align(nucl_num)

END_PROVIDER


BEGIN_PROVIDER  [ real, nucl_charge, (nucl_num) ]
  implicit none
  BEGIN_DOC
! Nuclear charge
  END_DOC

  nucl_charge = -1.0
  call get_nuclei_nucl_charge(nucl_charge)

  integer                        :: i
  do i=1,nucl_num
    if (nucl_charge(i) < 0.) then
      call abrt(irp_here,'Nuclear charges should be > 0')
    endif
  enddo
END_PROVIDER


BEGIN_PROVIDER [ real, nucl_coord,  (nucl_num_8,3) ]
  implicit none
  BEGIN_DOC
! Nuclear coordinates
  END_DOC

  nucl_coord = 0.
  ! The TREXIO (atom-major) layout only applies to non-worker processes.
  ! get_nuclei_nucl_coord() ignores use_trexio when is_worker is true: the
  ! dataserver sends the EZFIO layout (component-major, as stored in
  ! nuclei/nucl_coord). Keying this branch on use_trexio alone therefore
  ! scrambles the geometry of every worker as soon as use_trexio is true
  ! (E_nucl = 7.6664 instead of 9.1950 for water; E_loc = +48 instead of -75.6).
  if (use_trexio .and. (.not.is_worker)) then
    real, allocatable              :: buffer2(:,:)
    ! TREXIO stores /nucleus/nucleus_coord as an HDF5 dataset of dims
    ! (nucleus_num,3) read with a non-reversed memory dataspace
    ! (src/templates_hdf5/helpers_hdf5.c: H5Screate_simple(rank,size_read)),
    ! i.e. in C row-major order: x1,y1,z1,x2,y2,z2,...
    ! A Fortran (3,nucl_num) view has exactly that memory layout, so for a
    ! TREXIO read the buffer must be declared (3,nucl_num) and transposed.
    allocate (buffer2(3,nucl_num))
    call get_nuclei_nucl_coord(buffer2)
    do i=1,3
      do j=1,nucl_num
        nucl_coord(j,i) = buffer2(i,j)
      enddo
    enddo
    deallocate(buffer2)
  else
    real, allocatable              :: buffer(:,:)
    allocate (buffer(nucl_num,3))
    buffer = 0.
    call get_nuclei_nucl_coord(buffer)
    do i=1,3
      do j=1,nucl_num
        nucl_coord(j,i) = buffer(j,i)
      enddo
    enddo
    deallocate(buffer)
  end if
  integer                        :: i,j


END_PROVIDER


BEGIN_PROVIDER [ real, nucl_coord_transp, (8,nucl_num)
  implicit none
  BEGIN_DOC
! Transposed array of nucl_coord
  END_DOC
  integer                        :: i, k
  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    nucl_coord_transp = 0.
  endif

  !DIR$ VECTOR ALIGNED
  do i=1,nucl_num
    nucl_coord_transp(1,i) = nucl_coord(i,1)
    nucl_coord_transp(2,i) = nucl_coord(i,2)
    nucl_coord_transp(3,i) = nucl_coord(i,3)
  enddo
END_PROVIDER


 BEGIN_PROVIDER [ real, nucl_dist, (nucl_num_8,nucl_num) ]
&BEGIN_PROVIDER [ real, nucl_dist_vec_x, (nucl_num_8,nucl_num) ]
&BEGIN_PROVIDER [ real, nucl_dist_vec_y, (nucl_num_8,nucl_num) ]
&BEGIN_PROVIDER [ real, nucl_dist_vec_z, (nucl_num_8,nucl_num) ]
  implicit none
  BEGIN_DOC
! nucl_dist     : Nucleus-nucleus distances : nucl_dist(i,j) = |R_i-R_j|

! nucl_dist_vec : Nucleus-nucleus distances vectors
  END_DOC

  integer                        :: ie1, ie2, l
  integer,save                   :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    nucl_dist = 0.d0
    nucl_dist_vec_x = 0.d0
    nucl_dist_vec_y = 0.d0
    nucl_dist_vec_z = 0.d0
  endif

  do ie2 = 1,nucl_num
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do ie1 = 1,nucl_num
      nucl_dist_vec_x(ie1,ie2) = nucl_coord(ie1,1) - nucl_coord(ie2,1)
      nucl_dist_vec_y(ie1,ie2) = nucl_coord(ie1,2) - nucl_coord(ie2,2)
      nucl_dist_vec_z(ie1,ie2) = nucl_coord(ie1,3) - nucl_coord(ie2,3)
    enddo
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT (100)
    do ie1 = 1,nucl_num
      nucl_dist  (ie1,ie2) = nucl_dist_vec_x(ie1,ie2)*nucl_dist_vec_x(ie1,ie2) +&
          nucl_dist_vec_y(ie1,ie2)*nucl_dist_vec_y(ie1,ie2) +       &
          nucl_dist_vec_z(ie1,ie2)*nucl_dist_vec_z(ie1,ie2)
      nucl_dist(ie1,ie2) = sqrt(nucl_dist  (ie1,ie2))
      ASSERT (nucl_dist(ie1,ie2) > 0.)
    enddo
  enddo

END_PROVIDER

BEGIN_PROVIDER [ real, nucl_fitcusp_radius, (nucl_num) ]
  implicit none
  BEGIN_DOC
! Distance threshold for the fit
  END_DOC
  real                           :: def(nucl_num), factor
  integer                        :: k
  real, parameter                :: a = 1.74891
  real, parameter                :: b = 0.126057

  if (.not. do_nucl_fitcusp) then
    nucl_fitcusp_radius = 0.d0
    return
  endif

  do k=1,nucl_num
    nucl_fitcusp_radius(k) = nucl_fitcusp_factor/(a*nucl_charge(k)+b)
  enddo

  ! Avoid dummy atoms
  do k=1,nucl_num
    if (nucl_charge(k) < 5.d-1) then
      nucl_fitcusp_radius(k) = 0.d0
    endif
  enddo

END_PROVIDER


