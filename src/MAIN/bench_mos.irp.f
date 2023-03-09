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
      do j=1, ao_num
        print 10, j, i, 0, mo_value(j,i) 
      enddo
    end do
    do i=1,elec_num
      do j=1,ao_num
        print 10, j, i, 1, mo_grad_x(j,i) 
      enddo
    end do
    do i=1,elec_num
      do j=1,ao_num
        print 10, j, i, 2, mo_grad_y(j,i)
      enddo
    end do
    do i=1,elec_num
      do j=1,ao_num
        print 10, j, i, 3, mo_grad_z(j,i) 
      enddo
    end do
    do i=1,elec_num
      do j=1,ao_num
        print 10, j, i, 4, mo_lapl(j,i) 
      enddo
    end do
end

subroutine step2(imax)
  implicit none
  integer, intent(in)            :: imax
  integer                        :: i, j
  do j=1,imax
    PROVIDE mo_value
    TOUCH elec_coord
  enddo
end
