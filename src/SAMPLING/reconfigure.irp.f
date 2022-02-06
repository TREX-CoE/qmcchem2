subroutine reconfigure(ipos,w)
  implicit none
  integer, intent(inout)         :: ipos(*)
  double precision, intent(in)   :: w(*)
  
  integer                        :: kptab(walk_num), kmtab(walk_num)
  double precision               :: wp(walk_num), wm(walk_num)
  double precision               :: tmp
  
  double precision               :: dwalk_num
  
  tmp = 0.d0
  do k=1,walk_num
    tmp = tmp + w(k)
  enddo
  dwalk_num = dble(walk_num)/tmp

  integer                        :: kp, km
  kp=0
  km=0

  double precision               :: accup, accum
  accup = 0.d0
  accum = 0.d0

  integer                        :: k
  do k=1,walk_num
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

  if(kp+km /= walk_num) then
    print *,  kp, km
    call abrt(irp_here,'pb in reconfiguration +/-')
  endif

  if(dabs(accup-accum) > 1.d-11) then
    print *,  accup, accum
    call abrt(irp_here,'pb in reconfiguration')
  endif
  
  double precision               :: qmc_ranf, rand
  double precision               :: rando(walk_num)
  rand = qmc_ranf()
  do k=1,walk_num
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

