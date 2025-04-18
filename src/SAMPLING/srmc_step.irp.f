! Providers of *_srmc_block_walk
!==============================
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

t = """
 BEGIN_PROVIDER [ $T, $X_srmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_srmc_block_walk_kahan $D2 ]
&BEGIN_PROVIDER [ $T, $X_2_srmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_2_srmc_block_walk_kahan $D2 ]
 implicit none
 BEGIN_DOC
! SRMC averages of $X. Computed in E_loc_srmc_block_walk
 END_DOC
 $X_srmc_block_walk = 0.d0
 $X_srmc_block_walk_kahan = 0.d0
 $X_2_srmc_block_walk = 0.d0
 $X_2_srmc_block_walk_kahan = 0.d0
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



 BEGIN_PROVIDER [ double precision, E_loc_srmc_block_walk       ]
&BEGIN_PROVIDER [ double precision, E_loc_2_srmc_block_walk     ]
&BEGIN_PROVIDER [ double precision, E_loc_srmc_block_walk_kahan, (3) ]
&BEGIN_PROVIDER [ double precision, E_loc_2_srmc_block_walk_kahan, (3) ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Properties averaged over the block using the SRMC method
 END_DOC

 real, allocatable :: elec_coord_tmp(:,:,:)
 integer :: mod_align
 double precision :: E_loc_save(4,walk_num_dmc_max)
 double precision :: E_loc_save_tmp(4,walk_num_dmc_max)
 integer :: trapped_walk_tmp(walk_num_dmc_max)
 double precision :: psi_value_save(walk_num)
 double precision :: psi_value_save_tmp(walk_num)
 double precision :: srmc_weight(walk_num)
 double precision, allocatable :: psi_grad_psi_inv_save(:,:,:)
 double precision, allocatable :: psi_grad_psi_inv_save_tmp(:,:,:)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: E_loc_save
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: E_loc_save_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: srmc_weight
 allocate ( psi_grad_psi_inv_save(elec_num_8,3,walk_num) ,          &
       psi_grad_psi_inv_save_tmp(elec_num_8,3,walk_num) ,            &
       elec_coord_tmp(mod_align(elec_num+1),3,walk_num) )
 psi_value_save = 0.d0
 psi_value_save_tmp = 0.d0
 srmc_weight = 1.d0

 double precision, external     :: qmc_ranf

! Initialization
 if (vmc_algo /= t_Brownian) then
   call abrt(irp_here,'SRMC should run with Brownian algorithm')
 endif

 integer :: k, i_walk, i_step

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   !DIR$ VECTOR ALIGNED
   $X_srmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_srmc_block_walk_kahan = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_srmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_srmc_block_walk_kahan = 0.d0
 endif
"""
for p in properties:
 print (t.replace("$X",p[1]))
END_SHELL

 logical                      :: loop
 integer*8                    :: cpu0, cpu1, cpu2, count_rate, count_max

 loop = .True.
 call system_clock(cpu0, count_rate, count_max)
 cpu2 = cpu0

 block_weight = 0.d0

 real, external                 :: accep_rate
 double precision               :: delta

 logical :: first_loop
 first_loop = .True.

 do while (loop)

  ! Every walker makes a step
  do i_walk=1,walk_num

    integer                        :: i,j,l
    if (.not.first_loop) then
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
        E_loc = E_loc_save(1,i_walk)
      enddo
      SOFT_TOUCH elec_coord psi_grad_psi_inv_x psi_grad_psi_inv_y psi_grad_psi_inv_z psi_value E_loc
    else
      do l=1,3
        do i=1,elec_num+1
          elec_coord(i,l) = elec_coord_full(i,l,i_walk)
        enddo
      enddo
      TOUCH elec_coord
      psi_value_save(i_walk) = psi_value
      E_loc_save(:,i_walk) = E_loc + pseudo_stabilization
   endif

   double precision               :: p,q
   real                           :: delta_x
   logical                        :: accepted
   call brownian_step_reg(p,q,accepted,delta_x)
   if (accepted) then
      trapped_walk(i_walk) = 0
   else
      trapped_walk(i_walk) += 1
   endif

   if ( (trapped_walk(i_walk) < trapped_walk_max).and. &
        (psi_value * psi_value_save(i_walk) >= 0.d0) ) then

      ! 2-step
      delta = ((E_loc+pseudo_stabilization)+E_loc_save(1,i_walk))*0.5d0

!     ! 3-step
!     delta = (5.d0 * (E_loc+pseudo_stabilization) + 8.d0 * E_loc_save(1,i_walk) - E_loc_save(2,i_walk))/12.d0

!     ! 4-step
!     delta = (9.d0*(E_loc+pseudo_stabilization)+19.d0*E_loc_save(1,i_walk)- &
!            5.d0*E_loc_save(2,i_walk)+E_loc_save(3,i_walk))/24.d0

     ! 5-step
!     delta = -((-251.d0*(E_loc+pseudo_stabilization))-646.d0*E_loc_save(1,i_walk)+264.d0*E_loc_save(2,i_walk)-&
!          106.d0*E_loc_save(3,i_walk)+19.d0*E_loc_save(4,i_walk))/720.d0


     delta = (delta - E_ref)

     srmc_weight(i_walk) = dexp(-dtime_step*delta*p)

   else
     srmc_weight(i_walk) = 0.d0
     trapped_walk(i_walk) = 0
   endif

    ! Trick to avoid holes in DMC PES.
!    if (dabs(delta/E_ref) * time_step_sq > p * 0.5d0 ) then
!      srmc_weight(i_walk) = 0.d0
!    endif

   elec_coord(elec_num+1,1) += p*time_step
   elec_coord(elec_num+1,2)  = (E_loc+pseudo_stabilization)
   elec_coord(elec_num+1,3)  = srmc_weight(i_walk)
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
!  if (accepted) then
      E_loc_save(4,i_walk) = E_loc_save(3,i_walk)
      E_loc_save(3,i_walk) = E_loc_save(2,i_walk)
      E_loc_save(2,i_walk) = E_loc_save(1,i_walk)
      E_loc_save(1,i_walk) = (E_loc+pseudo_stabilization)
!  endif

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
     if (calc_$X) then
   ! Kahan's summation algorithm to compute these sums reducing the rounding error:
   !  $X_srmc_block_walk    += $X * srmc_pop_weight_mult * srmc_weight(i_walk)
   !  $X_2_srmc_block_walk  += $X_2 * srmc_pop_weight_mult * srmc_weight(i_walk)
   ! see http://en.wikipedia.org/wiki/Kahan_summation_algorithm

      $X_srmc_block_walk_kahan($D2 3) = $X * srmc_pop_weight_mult * srmc_weight(i_walk) - $X_srmc_block_walk_kahan($D2 1)
      $X_srmc_block_walk_kahan($D2 2) = $X_srmc_block_walk $D1  + $X_srmc_block_walk_kahan($D2 3)
      $X_srmc_block_walk_kahan($D2 1) = ($X_srmc_block_walk_kahan($D2 2) - $X_srmc_block_walk $D1 ) &
          - $X_srmc_block_walk_kahan($D2 3)
      $X_srmc_block_walk $D1  =  $X_srmc_block_walk_kahan($D2 2)


      $X_2_srmc_block_walk_kahan($D2 3) = $X_2 * srmc_pop_weight_mult * srmc_weight(i_walk) - $X_2_srmc_block_walk_kahan($D2 1)
      $X_2_srmc_block_walk_kahan($D2 2) = $X_2_srmc_block_walk $D1 + $X_2_srmc_block_walk_kahan($D2 3)
      $X_2_srmc_block_walk_kahan($D2 1) = ($X_2_srmc_block_walk_kahan($D2 2) - $X_2_srmc_block_walk $D1 ) &
          - $X_2_srmc_block_walk_kahan($D2 3)
      $X_2_srmc_block_walk $D1 =  $X_2_srmc_block_walk_kahan($D2 2)
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

   block_weight += srmc_pop_weight_mult * srmc_weight(i_walk)

  enddo ! i_walk

  ! Compute the new weight of the population
  double precision :: sum_weight
  sum_weight = 0.d0
  do k=1,walk_num
    sum_weight += srmc_weight(k)
  enddo
!  E_ref = E_ref - log(sum_weight/real(walk_num)) * 1.d0/srmc_projection_time

  ! Move to the next projection step
  if (srmc_projection > 0) then
    srmc_projection_step = mod(srmc_projection_step,srmc_projection)+1
  else
    srmc_projection_step = 1
  endif

  ! Eventually, recompute the weight of the population
  if (srmc_projection_step == 1) then
    srmc_pop_weight_mult = 1.d0
    do k=1,srmc_projection
      srmc_pop_weight_mult *= srmc_pop_weight(k)
    enddo
  endif

  ! Remove contribution of the old value of the weight at the new
  ! projection step
   srmc_pop_weight_mult *= 1.d0/srmc_pop_weight(srmc_projection_step)

  srmc_pop_weight(srmc_projection_step) = sum_weight/dble(walk_num)

  ! Update the running population weight
  srmc_pop_weight_mult *= srmc_pop_weight(srmc_projection_step)


  if (do_print_dmc_data) then
    do k=1,walk_num
      if (qmc_ranf() < 0.001) then
       print *, '--'
       do i=1,elec_num
         print *, elec_coord_full(i,1:3,k)
       enddo
       print *, 'w=', srmc_weight(k)
      endif
    enddo
  endif

! Reconfiguration
  integer :: ipos(walk_num)

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
    E_loc_save_tmp(:,k) = E_loc_save(:,k)
    trapped_walk_tmp(k) = trapped_walk(k)
    ipos(k) = k
  enddo

  call reconfigure(ipos,srmc_weight)

  integer :: ipm
  do k=1,walk_num
   ipm = ipos(k)
   do l=1,3
    do i=1,elec_num+1
     elec_coord_full(i,l,k) = elec_coord_tmp(i,l,ipm)
    enddo
    do i=1,elec_num
      psi_grad_psi_inv_save(i,l,k) = psi_grad_psi_inv_save_tmp(i,l,ipm)
    enddo
   enddo
   psi_value_save(k) = psi_value_save_tmp(ipm)
   E_loc_save(:,k) = E_loc_save_tmp(:,ipm)
   trapped_walk(k) = trapped_walk_tmp(ipm)
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

  SOFT_TOUCH elec_coord_full srmc_pop_weight_mult psi_value psi_grad_psi_inv_x psi_grad_psi_inv_y psi_grad_psi_inv_z elec_coord E_ref

  first_loop = .False.

 enddo ! while loop

 double precision :: factor
 factor = 1.d0/block_weight
 SOFT_TOUCH block_weight

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   $X_srmc_block_walk   *= factor
   $X_2_srmc_block_walk *= factor
 endif
"""
for p in properties:
 print (t.replace("$X",p[1]))
END_SHELL

 deallocate(elec_coord_tmp, psi_grad_psi_inv_save, psi_grad_psi_inv_save_tmp)

END_PROVIDER


BEGIN_PROVIDER [ double precision, srmc_pop_weight_mult ]
 implicit none
 BEGIN_DOC
! Population weight of SRMC
 END_DOC
 srmc_pop_weight_mult = srmc_pop_weight(srmc_projection)
END_PROVIDER

BEGIN_PROVIDER [ real, srmc_projection_time ]
 implicit none
 BEGIN_DOC
! SRMC project time in au
 END_DOC
 srmc_projection_time = 1.
 call get_simulation_srmc_projection_time(srmc_projection_time)
END_PROVIDER



 BEGIN_PROVIDER [ integer, srmc_projection ]
&BEGIN_PROVIDER [ integer, srmc_projection_step ]
 implicit none
 BEGIN_DOC
! Number of projection steps for SRMC
 END_DOC
 srmc_projection = int( srmc_projection_time/time_step)
 srmc_projection_step = 0
END_PROVIDER

BEGIN_PROVIDER [ double precision, srmc_pop_weight, (0:srmc_projection+1) ]
 implicit none
 BEGIN_DOC
! Population weight of SRMC
 END_DOC
 srmc_pop_weight(:) = 1.d0
END_PROVIDER


BEGIN_PROVIDER [ logical, do_print_dmc_data ]
 implicit none
 BEGIN_DOC
 ! If true, print in stdout the data to fit a Jastrow
 END_DOC
 do_print_dmc_data = .False.
END_PROVIDER

