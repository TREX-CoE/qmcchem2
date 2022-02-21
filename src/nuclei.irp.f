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

  nucl_charge = -1.d0
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
  real, allocatable              :: buffer(:,:)
  if (use_trexio) then
    allocate (buffer(3,nucl_num))
    call get_nucl_coord_trexio(buffer)
    do i=1,3
      do j=1,nucl_num
        nucl_coord(j,i) = buffer(i,j)
      enddo
    enddo
  else
    allocate (buffer(nucl_num,3))
    buffer = 0.
    call get_nuclei_nucl_coord(buffer)
    do i=1,3
      do j=1,nucl_num
        nucl_coord(j,i) = buffer(j,i)
      enddo
    enddo
  end if
  integer                        :: i,j

  deallocate(buffer)

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
    nucl_dist = 0.
    nucl_dist_vec_x = 0.
    nucl_dist_vec_y = 0.
    nucl_dist_vec_z = 0.
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
      nucl_fitcusp_radius(k) = 0.
    endif
  enddo

END_PROVIDER


