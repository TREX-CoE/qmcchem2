! Providers of *_fkmc_block_walk
!==============================
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

t = """
 BEGIN_PROVIDER [ $T, $X_fkmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_fkmc_block_walk_kahan $D2 ]
&BEGIN_PROVIDER [ $T, $X_2_fkmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_2_fkmc_block_walk_kahan $D2 ]
 implicit none
 BEGIN_DOC  
! fkMC averages of $X. Computed in E_loc_fkmc_block_walk
 END_DOC
 $X_fkmc_block_walk = 0.d0
 $X_fkmc_block_walk_kahan = 0.d0
 $X_2_fkmc_block_walk = 0.d0
 $X_2_fkmc_block_walk_kahan = 0.d0
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
    print(t.replace("$X",p[1]).replace("$T",p[0]).replace("$D1",D1).replace("$D2",D2))

END_SHELL



 BEGIN_PROVIDER [ double precision, E_loc_fkmc_block_walk       ]
&BEGIN_PROVIDER [ double precision, E_loc_2_fkmc_block_walk     ]
&BEGIN_PROVIDER [ double precision, E_loc_fkmc_block_walk_kahan, (3) ]
&BEGIN_PROVIDER [ double precision, E_loc_2_fkmc_block_walk_kahan, (3) ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Properties averaged over the block using the FKMC method
 END_DOC

  integer, parameter :: BIRTH=1, DEATH=2
  real, allocatable :: elec_coord_tmp(:,:,:)
  integer :: mod_align
  double precision :: E_loc_save(walk_num_dmc_max)
  double precision :: E_loc_save_tmp(walk_num_dmc_max)
  double precision :: psi_value_save(walk_num)
  double precision :: psi_value_save_tmp(walk_num)
  double precision :: fkmc_weight(walk_num)
  double precision :: delta(walk_num)
  double precision, allocatable :: psi_grad_psi_inv_save(:,:,:)
  double precision, allocatable :: psi_grad_psi_inv_save_tmp(:,:,:)
  double precision, allocatable :: fkmc_clock_tmp(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save_tmp
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: E_loc_save
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: E_loc_save_tmp
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save_tmp
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: fkmc_weight
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: fkmc_clock_tmp
  allocate ( psi_grad_psi_inv_save(elec_num_8,3,walk_num),           &
       psi_grad_psi_inv_save_tmp(elec_num_8,3,walk_num),             &
       elec_coord_tmp(mod_align(elec_num+1),3,walk_num),             &
       fkmc_clock_tmp(2,walk_num) )
  psi_value_save = 0.d0
  psi_value_save_tmp = 0.d0
  fkmc_weight = 1.d0

! Initialization
 if (vmc_algo /= t_Brownian) then
   call abrt(irp_here,'FKMC should run with Brownian algorithm')
 endif

 integer :: k, i_walk, i_step

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   !DIR$ VECTOR ALIGNED
   $X_fkmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_fkmc_block_walk_kahan = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_fkmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_fkmc_block_walk_kahan = 0.d0
 endif
"""
for p in properties:
 print (t.replace("$X",p[1]))
END_SHELL

 logical                        :: loop
 integer*8                      :: cpu0, cpu1, cpu2, count_rate, count_max

 loop = .True.
 call system_clock(cpu0, count_rate, count_max)
 cpu2 = cpu0

 block_weight = 0.d0

 real, external                 :: accep_rate
 double precision               :: thr

 thr = 2.d0/time_step_sq

 logical :: first_loop
 first_loop = .True.

 do while (loop)

  ! Every walker makes a step
  do i_walk=1,walk_num
    
    if (.not.first_loop) then
      integer                        :: i,j,l
      do l=1,3
        do i=1,elec_num+1
          elec_coord(i,l) = elec_coord_full(i,l,i_walk)
        enddo
        do i=1,elec_num
          psi_grad_psi_inv_x(i) = psi_grad_psi_inv_save(i,1,i_walk)
          psi_grad_psi_inv_y(i) = psi_grad_psi_inv_save(i,2,i_walk)
          psi_grad_psi_inv_z(i) = psi_grad_psi_inv_save(i,3,i_walk)
        enddo
        psi_value = psi_value_save(i_walk)
        E_loc = E_loc_save(i_walk)
      enddo
      SOFT_TOUCH elec_coord psi_grad_psi_inv_x psi_grad_psi_inv_y psi_grad_psi_inv_z psi_value E_loc
    else
      do l=1,3
        do i=1,elec_num+1
          elec_coord(i,l) = elec_coord_full(i,l,i_walk)
        enddo
      enddo
      TOUCH elec_coord
      E_loc_save(i_walk) = E_loc 
      psi_value_save(i_walk) = psi_value
    endif

   double precision               :: p,q
   real                           :: delta_x
   logical                        :: accepted
   call brownian_step(p,q,accepted,delta_x)

   if ( psi_value * psi_value_save(i_walk) >= 0.d0 ) then
     delta(i_walk) = ((E_loc+E_loc_save(i_walk))*0.5d0 - E_ref) * p
     if ( delta(i_walk) > thr ) then
       delta(i_walk) = thr
     else if ( delta(i_walk) < -thr ) then
       delta(i_walk) = -thr
     endif
     fkmc_weight(i_walk) = dexp(-dtime_step*delta(i_walk))
     elec_coord(elec_num+1,1) += p*time_step
     elec_coord(elec_num+1,2)  = E_loc
     elec_coord(elec_num+1,3)  = fkmc_weight(i_walk) 
     do l=1,3
        do i=1,elec_num+1
          elec_coord_full(i,l,i_walk) = elec_coord(i,l)
        enddo
     enddo
     do i=1,elec_num
       psi_grad_psi_inv_save(i,1,i_walk) = psi_grad_psi_inv_x(i)
       psi_grad_psi_inv_save(i,2,i_walk) = psi_grad_psi_inv_y(i)
       psi_grad_psi_inv_save(i,3,i_walk) = psi_grad_psi_inv_z(i)
     enddo

     psi_value_save(i_walk) = psi_value
     E_loc_save(i_walk) = E_loc

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
     if (calc_$X) then
   ! Kahan's summation algorithm to compute these sums reducing the rounding error:
   !  $X_fkmc_block_walk    += $X * fkmc_weight(i_walk)
   !  $X_2_fkmc_block_walk  += $X_2 * fkmc_weight(i_walk)
   ! see http://en.wikipedia.org/wiki/Kahan_summation_algorithm
   
      $X_fkmc_block_walk_kahan($D2 3) = $X * fkmc_weight(i_walk) - $X_fkmc_block_walk_kahan($D2 1)
      $X_fkmc_block_walk_kahan($D2 2) = $X_fkmc_block_walk $D1  + $X_fkmc_block_walk_kahan($D2 3)
      $X_fkmc_block_walk_kahan($D2 1) = ($X_fkmc_block_walk_kahan($D2 2) - $X_fkmc_block_walk $D1 ) &
          - $X_fkmc_block_walk_kahan($D2 3)
      $X_fkmc_block_walk $D1  =  $X_fkmc_block_walk_kahan($D2 2) 
   
   
      $X_2_fkmc_block_walk_kahan($D2 3) = $X_2 * fkmc_weight(i_walk) - $X_2_fkmc_block_walk_kahan($D2 1)
      $X_2_fkmc_block_walk_kahan($D2 2) = $X_2_fkmc_block_walk $D1 + $X_2_fkmc_block_walk_kahan($D2 3)
      $X_2_fkmc_block_walk_kahan($D2 1) = ($X_2_fkmc_block_walk_kahan($D2 2) - $X_2_fkmc_block_walk $D1 ) &
          - $X_2_fkmc_block_walk_kahan($D2 3)
      $X_2_fkmc_block_walk $D1 =  $X_2_fkmc_block_walk_kahan($D2 2) 
     endif
"""
for p in properties:
  if p[2] == "":
   D1 = ""
   D2 = ""
  else:
   D1 = "("+":"*(p[2].count(',')+1)+")"
   D2 = ":"*(p[2].count(',')+1)+","
  print(t.replace("$X",p[1]).replace("$D1",D1).replace("$D2",D2))

END_SHELL

    block_weight += fkmc_weight(i_walk)

   else
     fkmc_weight(i_walk) = 0.d0
     delta(i_walk) = 1.d5
   endif

  enddo

  ! Compute the new weight of the population
  double precision :: sum_weight
  sum_weight = 0.d0
  do k=1,walk_num
    sum_weight += fkmc_weight(k)
  enddo

  do k=1,walk_num
    do l=1,3
     do i=1,elec_num+1
      elec_coord_tmp(i,l,k) = elec_coord_full(i,l,k)
     enddo
     do i=1,elec_num
      psi_grad_psi_inv_save_tmp(i,l,k) = psi_grad_psi_inv_save(i,l,k)
     enddo
    enddo
    psi_value_save_tmp(k) = psi_value_save(k)
    E_loc_save_tmp(k) = E_loc_save(k)
    if (fkmc_weight(k) == 0.d0) then
      fkmc_clock(DEATH,k) = -1.d0
    endif
    if ( delta(k) <= 0.d0 ) then
      fkmc_clock_tmp(BIRTH,k) = fkmc_clock(BIRTH,k) +time_step * delta(k)
      fkmc_clock_tmp(DEATH,k) = fkmc_clock(DEATH,k)
    else
      fkmc_clock_tmp(BIRTH,k) = fkmc_clock(BIRTH,k)
      fkmc_clock_tmp(DEATH,k) = fkmc_clock(DEATH,k) -time_step * delta(k)
    endif
  enddo

! Reconfiguration
! ===============

  ! Identify first which walkers will be killed to place branched walkers there
  ! later
  
  double precision, external     :: qmc_ranf
  integer                        :: ipm, m
  integer                        :: killed(walk_num)

  m=1
  do k=1,walk_num
    fkmc_clock(DEATH,k) = fkmc_clock_tmp(DEATH,k) 
    if (fkmc_clock_tmp(DEATH,k) <= 0.d0) then
      killed(m) = k
      m += 1
      fkmc_clock(DEATH,k) = -dlog(qmc_ranf())
      fkmc_clock(BIRTH,k) = -dlog(qmc_ranf())
      ipm = k
      do while (ipm == k)
        ipm = 1 + int (walk_num*qmc_ranf())
      enddo
      do l=1,3
        do i=1,elec_num+1
        elec_coord_full(i,l,k) = elec_coord_tmp(i,l,ipm)
        enddo
        do i=1,elec_num
          psi_grad_psi_inv_save(i,l,k) = psi_grad_psi_inv_save_tmp(i,l,ipm)
        enddo
      enddo
      psi_value_save(k) = psi_value_save_tmp(ipm)
      E_loc_save(k) = E_loc_save_tmp(ipm)
    endif
  enddo
  killed(m) = 0

  m=1
  do k=1,walk_num
    fkmc_clock(BIRTH,k) = fkmc_clock_tmp(BIRTH,k) 
    if (fkmc_clock_tmp(BIRTH,k) <= 0.d0) then
      fkmc_clock(BIRTH,k) = -dlog(qmc_ranf())
      if (killed(m) == 0) then
        ipm = k
        do while (ipm == k)
          ipm = 1 + int (walk_num*qmc_ranf())
        enddo
      else
        ipm = killed(m)
        m +=1
      endif
      fkmc_clock(BIRTH,ipm) = -dlog(qmc_ranf())
      fkmc_clock(DEATH,ipm) = -dlog(qmc_ranf())
      do l=1,3
        do i=1,elec_num+1
        elec_coord_full(i,l,ipm) = elec_coord_tmp(i,l,k)
        enddo
        do i=1,elec_num
          psi_grad_psi_inv_save(i,l,ipm) = psi_grad_psi_inv_save_tmp(i,l,k)
        enddo
      enddo
      psi_value_save(ipm) = psi_value_save_tmp(k)
      E_loc_save(ipm) = E_loc_save_tmp(k)
      
    endif

  enddo

  

  call system_clock(cpu1, count_rate, count_max)
  if (cpu1 < cpu0) then
    cpu1 = cpu1+cpu0
  endif
  loop = dble(cpu1-cpu0)/dble(count_rate) < block_time
  if (cpu1-cpu2 > count_rate) then
    integer                        :: do_run
    call get_running(do_run)
    loop = loop.and.(do_run == t_Running)
    cpu2 = cpu1
  endif

! Update E_ref to take into account the weight of the population
  E_ref -= dlog(sum_weight / dble(walk_num) ) / time_step 
  SOFT_TOUCH elec_coord_full E_ref

  first_loop = .False.

 enddo

 double precision :: factor
 factor = 1.d0/block_weight
 SOFT_TOUCH block_weight

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   $X_fkmc_block_walk   *= factor
   $X_2_fkmc_block_walk *= factor
 endif
"""
for p in properties:
 print (t.replace("$X",p[1]))
END_SHELL

 deallocate ( elec_coord_tmp, psi_grad_psi_inv_save, psi_grad_psi_inv_save_tmp, &
   fkmc_clock_tmp )

END_PROVIDER



BEGIN_PROVIDER [ double precision, fkmc_clock, (2,walk_num) ]
 implicit none
 BEGIN_DOC
 ! Branching clocks for the FKMC algotithm. (1,:) is the birth clock and
 ! (2,:) is the death clock.
 END_DOC
 integer :: i
 double precision, external :: qmc_ranf
 do i=1, walk_num
   fkmc_clock(1,i) = -dlog(qmc_ranf())
   fkmc_clock(2,i) = -dlog(qmc_ranf())
 enddo

END_PROVIDER

