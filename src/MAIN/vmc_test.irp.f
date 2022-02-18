program vmc_test
  real :: t1,t2
  print *,  'Ndet=',det_num
  print *,  'Ndet alpha beta =',det_alpha_num, det_beta_num
  if (do_prepare) then
    stop 'No walkers'
  endif
  print *,  'E_loc = ', E_loc
  call step2
  call ezfio_finish
end

subroutine step2
  implicit none
  real :: accep_rate
  print *,  '---'
  print *,  '<E_loc> = ', E_loc_block_walk
  print *,  '<E_loc_2> = ', E_loc_2_block_walk
  print *,  'w = ', block_weight
  print *,  'Accept', accep_rate()
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


