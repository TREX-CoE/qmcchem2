 use qmckl
! Providers of *_vmc_block_walk
!==============================
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

t = """
 BEGIN_PROVIDER [ $T, $X_vmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_vmc_block_walk_kahan $D2 ]
&BEGIN_PROVIDER [ $T, $X_2_vmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_2_vmc_block_walk_kahan $D2 ]
 implicit none
 BEGIN_DOC
! VMC averages of $X
 END_DOC
END_PROVIDER
"""
for p in properties:
  if p[1] != 'e_loc':
    if p[2] == "":
      D1 = ""
      D2 = ", (3)"
    else:
      D1 = ", ("+p[2][1:-1]+")"
      D2 = ", ("+p[2][1:-1]+",3)"
    print (t.replace("$X",p[1]).replace("$T",p[0]).replace("$D1",D1).replace("$D2",D2))
END_SHELL

 BEGIN_PROVIDER [ double precision, E_loc_vmc_block_walk ]
&BEGIN_PROVIDER [ double precision, E_loc_2_vmc_block_walk ]
&BEGIN_PROVIDER [ double precision, E_loc_vmc_block_walk_kahan, (3) ]
&BEGIN_PROVIDER [ double precision, E_loc_2_vmc_block_walk_kahan, (3) ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Properties averaged over the block per walker using the VMC method
 END_DOC

 integer :: i_walk

 PROVIDE time_step

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
  if (calc_$X) then
   !DIR$ VECTOR ALIGNED
   $X_vmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_vmc_block_walk_kahan = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_vmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_vmc_block_walk_kahan = 0.d0
   $X_min = huge(1.)
   $X_max =-huge(1.)
  endif
"""
for p in properties:
 print (t.replace("$X",p[1]))
END_SHELL

 double precision :: dnorm
 block_weight = 0.d0
 do i_walk=1,walk_num
   integer :: i,j,l
   if (i_walk > 1) then
     do l=1,3
      do i=1,elec_num+1
       elec_coord(i,l) = elec_coord_full(i,l,i_walk)
      enddo
     enddo
     call update_qmckl_coord()
     TOUCH elec_coord
   endif

   logical                        :: loop
   integer*8                      :: cpu0, cpu1, cpu2, count_rate, count_max
   loop = .True.

   call system_clock(cpu0, count_rate, count_max)
   cpu2 = cpu0
   do while (loop)
     double precision               :: p,q
     real                           :: delta_x
     logical                        :: accepted
     if (vmc_algo == t_Brownian) then
       call brownian_step(p,q,accepted,delta_x)
     else if (vmc_algo == t_Langevin) then
       call langevin_step(p,q,accepted,delta_x)
     endif
     elec_coord(elec_num+1,1) += p*time_step
     elec_coord(elec_num+1,2) = E_loc
     elec_coord(elec_num+1,3) = 1.
     block_weight += 1.d0
   

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
     if (calc_$X) then
   ! Kahan's summation algorithm to compute these sums reducing the rounding error:
   !  $X_vmc_block_walk $D1   += $X
   !  $X_2_vmc_block_walk $D1 += $X_2
   ! see http://en.wikipedia.org/wiki/Kahan_summation_algorithm

      $X_vmc_block_walk_kahan($D2 3) = $X - $X_vmc_block_walk_kahan($D2 1)
      $X_vmc_block_walk_kahan($D2 2) = $X_vmc_block_walk $D1  + $X_vmc_block_walk_kahan($D2 3)
      $X_vmc_block_walk_kahan($D2 1) = ($X_vmc_block_walk_kahan($D2 2) - $X_vmc_block_walk $D1 ) &
          - $X_vmc_block_walk_kahan($D2 3)
      $X_vmc_block_walk $D1  =  $X_vmc_block_walk_kahan($D2 2)


      $X_2_vmc_block_walk_kahan($D2 3) = $X_2 - $X_2_vmc_block_walk_kahan($D2 1)
      $X_2_vmc_block_walk_kahan($D2 2) = $X_2_vmc_block_walk $D1 + $X_2_vmc_block_walk_kahan($D2 3)
      $X_2_vmc_block_walk_kahan($D2 1) = ($X_2_vmc_block_walk_kahan($D2 2) - $X_2_vmc_block_walk $D1 ) &
          - $X_2_vmc_block_walk_kahan($D2 3)
      $X_2_vmc_block_walk $D1 =  $X_2_vmc_block_walk_kahan($D2 2)
     endif
"""
for p in properties:
  if p[2] == "":
   D1 = ""
   D2 = ""
  else:
   D1 = "("+":"*(p[2].count(',')+1)+")"
   D2 = ":"*(p[2].count(',')+1)+","
  print (t.replace("$X",p[1]).replace("$D1",D1).replace("$D2",D2))

END_SHELL


      call system_clock(cpu1, count_rate, count_max)
      if (cpu1 < cpu0) then
        cpu1 = cpu1+cpu0
      endif
      loop = dble(cpu1-cpu0)*dble(walk_num)/dble(count_rate) < block_time
      if (cpu1-cpu2 > count_rate) then
        integer                        :: do_run
        call get_running(do_run)
        loop = loop.and.(do_run == t_Running)
        cpu2 = cpu1
      endif

   enddo ! while (loop)

   do l=1,3
    do i=1,elec_num+1
     elec_coord_full(i,l,i_walk) = elec_coord(i,l)
    enddo
   enddo
   
 enddo

 double precision               :: factor
 factor = 1.d0/block_weight
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   $X_vmc_block_walk *= factor
   $X_2_vmc_block_walk *= factor
 endif
"""
for p in properties:
 print  (t.replace("$X",p[1]))
END_SHELL

 SOFT_TOUCH elec_coord_full block_weight

END_PROVIDER


