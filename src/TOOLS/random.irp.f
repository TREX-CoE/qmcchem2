double precision function qmc_ranf()
! L'Ecuyer, P. (1999) `Tables of maximally equidistributed combined LFSR
! generators', Math. of Comput., 68, 261-269.
 implicit none
 integer*8 :: b(2)
 b(1) = SHIFTR( IEOR( ISHFT(seed(1),1), seed(1)), 53)
 b(2) = SHIFTL( IAND(seed(1),-2_8), 10)
 seed(1) = IEOR( b(2), b(1))

 b(1) = SHIFTR( IEOR( ISHFT(seed(2),24), seed(2)), 50)
 b(2) = SHIFTL( IAND(seed(2),-512_8), 5)
 seed(2) = IEOR( b(2), b(1))

 b(1) = SHIFTR( IEOR( ISHFT(seed(3),3), seed(3)), 23)
 b(2) = SHIFTL( IAND(seed(3),-4096_8), 29)
 seed(3) = IEOR( b(2), b(1))

 b(1) = SHIFTR( IEOR( ISHFT(seed(4),5), seed(4)), 24)
 b(2) = SHIFTL( IAND(seed(4),-131072_8), 23)
 seed(4) = IEOR( b(2), b(1))

 b(1) = SHIFTR( IEOR( ISHFT(seed(5),3), seed(5)), 33)
 b(2) = SHIFTL( IAND(seed(5),-8388608_8), 8)
 seed(5) = IEOR( b(2), b(1))

 qmc_ranf = IEOR( IEOR( IEOR( IEOR(seed(1),seed(2)), seed(3)), &
   seed(4)), seed(5)) * 5.4210108624275221D-20 + 0.5D0
 ASSERT ( qmc_ranf >= 0.d0 )
 ASSERT ( qmc_ranf <= 1.d0 )

end

subroutine ranf_array(isize,res)
 implicit none
 integer :: isize
 double precision :: res(isize)
 integer :: i
 double precision :: qmc_ranf

 do i=1,isize
   res(i) = qmc_ranf()
 enddo
end

BEGIN_PROVIDER [ logical, deterministic ]
 implicit none
 BEGIN_DOC
 ! Only false in qmc program
 END_DOC
 deterministic = .True.
END_PROVIDER

BEGIN_PROVIDER  [ integer*8, seed, (5) ]
  implicit none
  BEGIN_DOC  
! Seeds data
! Initialized by init_random
  END_DOC
  integer                        :: iargc
  integer*8                      :: i,j
  integer*4                      :: clock(33)
  double precision               :: r
  integer*8                      :: pid8

  if (deterministic) then
    do i=1,size(clock)
      clock(i) = i
    enddo
    call random_seed(put=clock)

  else
    read(current_PID,*) pid8
    pid8 = iand( shiftl(pid8, 32), pid8)
    do i=1,size(clock)
      clock(i) = i
    enddo
    call system_clock(count=clock(1))
    call random_seed(put=clock)
  endif

    do i=1,5
      call random_number(r)
      seed(i) = (r-0.5d0)*huge(1_8)
      seed(i) = ieor( seed(i), pid8)
      do j=1,16
        seed(i) = shiftl(seed(i),1)+1
      enddo
    enddo

END_PROVIDER

subroutine gauss_array(isize,res)
  implicit none
  include '../constants.F'
  integer isize
  double precision res(isize)

  double precision u1(isize),u2(isize)
  integer i

  call ranf_array(isize,u1)
  call ranf_array(isize,u2)
  do i=1,isize
    res(i)=sqrt(-2.d0*log(u1(i)))*cos(dtwo_pi*u2(i))
  enddo
end

double precision function gauss()
  implicit none
  include '../constants.F'
  double precision :: qmc_ranf
  double precision :: u1,u2
  u1=qmc_ranf()
  u2=qmc_ranf()
  gauss=dsqrt(-2.d0*dlog(u1))*dcos(dfour_pi*u2)
end

