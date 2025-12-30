subroutine reconfigure_simple(ipos,w, isize)
  implicit none
  integer, intent(in)            :: isize
  integer, intent(inout)         :: ipos(isize)
  double precision, intent(in)   :: w(isize)

  integer                        :: i, j, k, n

  integer :: ipos_new(2*isize)
  double precision :: r, u(2*isize)
  double precision, external     :: qmc_ranf

  n = isize
  ipos_new(1:isize) = ipos(1:isize)
  do k=1,isize
    r = qmc_ranf()
    if (w(k) > 1.d0) then
      if ( r < w(k) - 1.d0 ) then
        n = n+1
        ipos_new(n) = k
      endif
    else
      if ( r > w(k) ) then
        ipos_new(k) = ipos_new(n)
        n = n-1
      endif
    endif
  enddo

  if (n == isize) then
    ipos(1:isize) = ipos_new(1:isize)
  else
    ! Randomly shuffle ipos_new
    call random_number(u(1:n))
    do k=1, n
     j = 1 + floor(dble(n)*u(k))
     i = ipos_new(k)
     ipos_new(k) = ipos_new(j)
     ipos_new(j) = i
    end do

    if (n > isize) then
      ipos(1:isize) = ipos_new(1:isize)
    else
      ipos(1:n) = ipos_new(1:n)
      ipos(n+1:isize) = ipos_new(1:isize-n)
    end if
  end if
end


subroutine reconfigure(ipos,w,isize)
  implicit none
  integer, intent(in)            :: isize
  integer, intent(inout)         :: ipos(isize)
  double precision, intent(in)   :: w(isize)

  integer                        :: i, j, k, n

  double precision, external     :: qmc_ranf
  integer                        :: kptab(isize), kmtab(isize)
  double precision               :: wp(isize), wm(isize)
  double precision               :: tmp

  double precision               :: dwalk_num

  tmp = 0.d0
  do k=1,isize
    tmp = tmp + w(k)
  enddo
  dwalk_num = dble(isize)/tmp

  integer                        :: kp, km
  kp=0
  km=0

  double precision               :: accup, accum
  accup = 0.d0
  accum = 0.d0

  do k=1,isize
    tmp = dwalk_num*w(k)-1.d0
    if (tmp >= 0.d0) then
      kp = kp+1
      wp(kp) = dabs(tmp)
      accup = accup + wp(kp)
      kptab(kp) = k
    else
      km = km+1
      wm(km) = dabs(tmp)
      accum = accum + wm(km)
      kmtab(km) = k
    endif
  enddo

  if(dabs(accup-accum) > 1.d-11) then
    print *,  accup, accum
    call abrt(irp_here,'pb in reconfiguration')
  endif

  double precision               :: rand
  double precision               :: rando(isize)
  rand = qmc_ranf()
  do k=1,isize
    rando(k) = dble(k-1)+rand
  enddo

  double precision               :: averageconf, current
  integer                        :: kcp
  integer                        :: kadd, kremove

  averageconf = accup
  kcp = 1
  rand = rando(kcp)

  do while (rand < averageconf)
    k=1
    current=wm(k)
    do while (rand > current)
      k = k+1
      current = current + wm(k)
    enddo
    kremove = kmtab(k)

    k=1
    current=wp(k)
    do while (rand > current)
      k = k+1
      current = current + wp(k)
    enddo
    kadd = kptab(k)

    ipos(kremove) = kadd
    kcp = kcp + 1
    rand = rando(kcp)
  enddo

end

