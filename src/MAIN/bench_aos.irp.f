program qmcchem_info
  implicit none
  PROVIDE ezfio_filename
  double precision               :: cpu0, cpu1
  character*(8)                  :: str_n
  integer                        :: iargc
  integer                        :: imax
  if (command_argument_count() > 1) then
    call get_command_argument(2,str_n)
    read(str_n,*) imax
  else
    imax = 100
  endif
  call step1
  print *,  '---'

  call wall_time (cpu0)
  call step2(imax)
  call wall_time (cpu1)
  print *,  'QMC=Chem : ', 1000.*(cpu1-cpu0)/float(imax)

  call wall_time (cpu0)
  call step3(imax)
  call wall_time (cpu1)
  print *,  'QMCkl    : ', 1000.*(cpu1-cpu0)/float(imax)
end

subroutine step1
  implicit none
  integer                        :: i, j
 10 format (3I4,2E20.10)
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1, ao_value_non_zero_idx(0)
        print 10, ao_value_non_zero_idx(j), i, 1, ao_value_block(j), &
          qmckl_ao_vgl(ao_value_non_zero_idx(j), 1, i)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1, ao_value_non_zero_idx(0)
        print 10, ao_value_non_zero_idx(j), i, 2, ao_grad_block_x(j), &
          qmckl_ao_vgl(ao_value_non_zero_idx(j), 2, i)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1, ao_value_non_zero_idx(0)
        print 10, ao_value_non_zero_idx(j), i, 3, ao_grad_block_y(j), &
          qmckl_ao_vgl(ao_value_non_zero_idx(j), 3, i)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1, ao_value_non_zero_idx(0)
        print 10, ao_value_non_zero_idx(j), i, 4, ao_grad_block_z(j), &
          qmckl_ao_vgl(ao_value_non_zero_idx(j), 4, i)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1, ao_value_non_zero_idx(0)
        print 10, ao_value_non_zero_idx(j), i, 5, ao_lapl_block(j), &
          qmckl_ao_vgl(ao_value_non_zero_idx(j), 5, i)
      enddo
    end do
end

subroutine step2(imax)
  implicit none
  integer, intent(in)            :: imax
  integer                        :: i, j
  do j=1,imax
    do i=1,elec_num
      if (i>1) then
        ao_elec = i
        TOUCH ao_elec
      endif
      PROVIDE ao_value_block
    end do
    TOUCH elec_coord_full
  enddo
end

subroutine step3(imax)
  implicit none
  integer, intent(in)            :: imax
  integer                        :: i, j
  do j=1,imax
    PROVIDE qmckl_ao_vgl
    TOUCH elec_coord
  enddo
end
