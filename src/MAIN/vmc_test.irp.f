program vmc_test
  real :: t1,t2
  print *,  'Ndet=',det_num
  print *,  'Ndet alpha beta =',det_alpha_num, det_beta_num
  if (do_prepare) then
    stop 'No walkers'
  endif
  print *, jast_value
  print *,  'E_loc = ', E_loc
  print *,  'Eloc_Jpsi = ', eloc_jpsi
!  print *, transpose(elec_coord(1:elec_num,1:3))
  call step2
  call ezfio_finish
end

subroutine step2
  implicit none
  real :: accep_rate
  integer, parameter :: NMAX=1000
  double precision :: E(NMAX), EJ(NMAX), ave, err, w
  integer :: i
  w = 0.d0
  do i=1,NMAX
    E(i) = E_loc_block_walk
    EJ(i) = Eloc_Jpsi_block_walk
    w += block_weight
    TOUCH elec_coord

    ave = sum(E(1:i))/dble(i)
    err = dsqrt(sum( (E(1:i)-ave)**2 ) / dble(i*i-i))

    print *, E_loc - Eloc_Jpsi
    print *,  '---'
    print *,  '<E_loc> = ', ave, '+/-', real(err)
    print *,  '<E_loc_2> = ', E_loc_2_block_walk
    ave = sum(EJ(1:i))/dble(i)
    err = dsqrt(sum( (EJ(1:i)-ave)**2 ) / dble(i*i-i))
    print *,  '<E_loc_Jpsi> = ', ave, '+/-', real(err)
    print *,  '<E_loc_Jpsi2> = ', Eloc_Jpsi_2_block_walk
    print *,  'w = ', w
    print *,  'Accept', accep_rate()
  end do
  print *,  ''
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

derivlist = []

for p in properties:
  t = """
    if (calc_$X) then
      PROVIDE $X_block_walk
      PROVIDE $X_2_block_walk
    endif
    """
  print(t.replace("$X",p[1]))
END_SHELL
end


