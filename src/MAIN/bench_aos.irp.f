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
  call cpu_time (cpu0)
  call step2(imax)
  call cpu_time (cpu1)
  print *,  'Time for the calculation of 1 step (ms)   : ', 1000.*(cpu1-cpu0)/float(imax)
end

subroutine step1
  implicit none
  integer                        :: i, j
 10 format (3I4,2E20.10)
!    print *, elec_coord(1:elec_num,1:3)
    print *, ao_nucl_sort_idx(:)
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1, ao_value_non_zero_idx(0)
        print 10, ao_value_non_zero_idx(j), i, 0, ao_value_block(j), elec_coord(i,1)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1,ao_num
        print 10, ao_value_non_zero_idx(j), i, 1, ao_grad_block_x(j)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1,ao_num
        print 10, ao_value_non_zero_idx(j), i, 2, ao_grad_block_y(j)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1,ao_num
        print 10, ao_value_non_zero_idx(j), i, 3, ao_grad_block_z(j)
      enddo
    end do
    do i=1,elec_num
      ao_elec = i
      TOUCH ao_elec
      do j=1,ao_num
        print 10, ao_value_non_zero_idx(j), i, 4, ao_lapl_block(j)
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
    TOUCH elec_coord
  enddo
end
