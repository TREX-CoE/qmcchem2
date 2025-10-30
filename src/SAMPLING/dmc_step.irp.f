! Providers of *_dmc_block_walk
!==============================
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

t = """
 BEGIN_PROVIDER [ $T, $X_dmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_dmc_block_walk_kahan $D2 ]
&BEGIN_PROVIDER [ $T, $X_2_dmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_2_dmc_block_walk_kahan $D2 ]
 implicit none
 BEGIN_DOC  
! DMC averages of $X. Computed in E_loc_dmc_block_walk
 END_DOC
 $X_dmc_block_walk = 0.d0
 $X_dmc_block_walk_kahan = 0.d0
 $X_2_dmc_block_walk = 0.d0
 $X_2_dmc_block_walk_kahan = 0.d0
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



 BEGIN_PROVIDER [ double precision, E_loc_dmc_block_walk       ]
&BEGIN_PROVIDER [ double precision, E_loc_2_dmc_block_walk     ]
&BEGIN_PROVIDER [ double precision, E_loc_dmc_block_walk_kahan, (3) ]
&BEGIN_PROVIDER [ double precision, E_loc_2_dmc_block_walk_kahan, (3) ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Properties averaged over the block using the DMC method
 END_DOC

 real, allocatable :: elec_coord_tmp(:,:,:)
 integer :: mod_align
 double precision :: E_loc_save(4,walk_num_dmc_max)
 double precision :: E_loc_save_tmp(4,walk_num_dmc_max)
 integer :: trapped_walk_tmp(walk_num_dmc_max)
 double precision :: psi_value_save(walk_num_dmc_max)
 double precision :: psi_value_save_tmp(walk_num_dmc_max)
 double precision :: dmc_weight(walk_num_dmc_max)
 double precision, allocatable :: psi_grad_psi_inv_save(:,:,:)
 double precision, allocatable :: psi_grad_psi_inv_save_tmp(:,:,:)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: E_loc_save
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: E_loc_save_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dmc_weight
 allocate ( psi_grad_psi_inv_save(elec_num_8,3,walk_num_dmc_max) ,          &
       psi_grad_psi_inv_save_tmp(elec_num_8,3,walk_num_dmc_max) ,            &
       elec_coord_tmp(mod_align(elec_num+1),3,walk_num_dmc_max) )
 psi_value_save = 0.d0
 psi_value_save_tmp = 0.d0
 dmc_weight = 1.d0

 double precision, external     :: qmc_ranf

! Initialization
 if (vmc_algo /= t_Brownian) then
   call abrt(irp_here,'DMC should run with Brownian algorithm')
 endif

 integer :: k, i_walk, i_step

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   !DIR$ VECTOR ALIGNED
   $X_dmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_dmc_block_walk_kahan = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_dmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_dmc_block_walk_kahan = 0.d0
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
 double precision               :: delta, E0

 logical :: first_loop
 first_loop = .True.

 E0 = E_ref

 do while (loop)

  ! Every walker makes a step
  do i_walk=1,walk_num_dmc

    integer                        :: i,j,l
    if (.not.first_loop) then
      do l=1,3
        do i=1,elec_num+1
          elec_coord(i,l) = elec_coord_full_dmc(i,l,i_walk)
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
          elec_coord(i,l) = elec_coord_full_dmc(i,l,i_walk)
        enddo
      enddo
      TOUCH elec_coord
      psi_value_save(i_walk) = psi_value
      E_loc_save(:,i_walk) = E_loc
   endif

   double precision               :: p,q
   real                           :: delta_x
   logical                        :: accepted
   call brownian_step(p,q,accepted,delta_x)
   if (accepted) then
      trapped_walk(i_walk) = 0
   else
      trapped_walk(i_walk) += 1
   endif

   if ( (trapped_walk(i_walk) < trapped_walk_max).and. &
        (psi_value * psi_value_save(i_walk) >= 0.d0) ) then

!     ! 2-step
!     delta = (E_loc+E_loc_save(1,i_walk))*0.5d0

     ! 4-step
     delta = (9.d0*E_loc+19.d0*E_loc_save(1,i_walk)- &
              5.d0*E_loc_save(2,i_walk)+E_loc_save(3,i_walk))/24.d0

     delta = (delta - E0)*p

     if (delta >= 0.d0) then
       dmc_weight(i_walk) = dexp(-dtime_step*delta)
     else
       dmc_weight(i_walk) = min(.05d0*walk_num,dexp(-dtime_step*delta))
!       dmc_weight(i_walk) = max(2.d0-dexp(dtime_step*delta), 0.d0)
     endif

   else
     dmc_weight(i_walk) = 0.d0
     trapped_walk(i_walk) = 0
   endif

    ! Trick to avoid holes in DMC PES.
!   if (dabs(delta/E_ref) * time_step_sq > p * 0.5d0 ) then
!     dmc_weight(i_walk) = 0.d0
!   endif

   elec_coord(elec_num+1,1) += p*time_step
   elec_coord(elec_num+1,2)  = E_loc
   elec_coord(elec_num+1,3)  = dmc_weight(i_walk)
   do l=1,3
      do i=1,elec_num+1
        elec_coord_full_dmc(i,l,i_walk) = elec_coord(i,l)
      enddo
   enddo
   do i=1,elec_num
     psi_grad_psi_inv_save(i,1,i_walk) = psi_grad_psi_inv_x(i)
     psi_grad_psi_inv_save(i,2,i_walk) = psi_grad_psi_inv_y(i)
     psi_grad_psi_inv_save(i,3,i_walk) = psi_grad_psi_inv_z(i)
   enddo

   psi_value_save(i_walk) = psi_value
   if (accepted) then
      E_loc_save(4,i_walk) = E_loc_save(3,i_walk)
      E_loc_save(3,i_walk) = E_loc_save(2,i_walk)
      E_loc_save(2,i_walk) = E_loc_save(1,i_walk)
      E_loc_save(1,i_walk) = E_loc
   endif

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
     if (calc_$X) then
   ! Kahan's summation algorithm to compute these sums reducing the rounding error:
   !  $X_dmc_block_walk    += $X * dmc_weight(i_walk)
   !  $X_2_dmc_block_walk  += $X_2 * dmc_weight(i_walk)
   ! see http://en.wikipedia.org/wiki/Kahan_summation_algorithm
   
      $X_dmc_block_walk_kahan($D2 3) = $X * dmc_weight(i_walk) - $X_dmc_block_walk_kahan($D2 1)
      $X_dmc_block_walk_kahan($D2 2) = $X_dmc_block_walk $D1  + $X_dmc_block_walk_kahan($D2 3)
      $X_dmc_block_walk_kahan($D2 1) = ($X_dmc_block_walk_kahan($D2 2) - $X_dmc_block_walk $D1 ) &
          - $X_dmc_block_walk_kahan($D2 3)
      $X_dmc_block_walk $D1  =  $X_dmc_block_walk_kahan($D2 2) 
   
   
      $X_2_dmc_block_walk_kahan($D2 3) = $X_2 * dmc_weight(i_walk) - $X_2_dmc_block_walk_kahan($D2 1)
      $X_2_dmc_block_walk_kahan($D2 2) = $X_2_dmc_block_walk $D1 + $X_2_dmc_block_walk_kahan($D2 3)
      $X_2_dmc_block_walk_kahan($D2 1) = ($X_2_dmc_block_walk_kahan($D2 2) - $X_2_dmc_block_walk $D1 ) &
          - $X_2_dmc_block_walk_kahan($D2 3)
      $X_2_dmc_block_walk $D1 =  $X_2_dmc_block_walk_kahan($D2 2) 
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

   block_weight += dmc_weight(i_walk)

  enddo ! i_walk

  ! Population control
  double precision :: sum_weight
  sum_weight = 0.d0
  do k=1,walk_num_dmc
    sum_weight += dmc_weight(k)
  enddo
  E0 = E_ref - log(sum_weight/dble(walk_num)) * 0.1d0 /dtime_step

! Branching
  integer                        :: walk_num_dmc_new
  integer                        :: ipos(walk_num_dmc_max)

  do k=1,walk_num_dmc
    do l=1,3
     do i=1,elec_num+1
      elec_coord_tmp(i,l,k) = elec_coord_full_dmc(i,l,k)
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

! call reconfigure(ipos,dmc_weight)
! walk_num_dmc = walk_num
  
  walk_num_dmc_new = walk_num_dmc
  double precision               :: r, u
  do k=1,walk_num_dmc
    u = dmc_weight(k)
    do while (u >= 1.d0)
      walk_num_dmc_new = walk_num_dmc_new+1
      if (walk_num_dmc_new > walk_num_dmc_max) exit
      ipos(walk_num_dmc_new) = k
      u = u-1.d0
    end do
    if (walk_num_dmc_new > walk_num_dmc_max) exit
    r = qmc_ranf()
    if ( r < u ) then
      ipos(k) = ipos(walk_num_dmc_new)
      walk_num_dmc_new = walk_num_dmc_new+1
    endif
    if (walk_num_dmc_new > walk_num_dmc_max) exit
  enddo
  if (walk_num_dmc_new == 0) then
    call abrt(irp_here,'Population disappeared')
  end if
  if (walk_num_dmc_new > walk_num_dmc_max) then
    call abrt(irp_here,'Population explosion')
  end if
  walk_num_dmc = walk_num_dmc_new


  integer :: ipm
  do k=1,walk_num_dmc
   ipm = ipos(k)
   do l=1,3
    do i=1,elec_num+1
     elec_coord_full_dmc(i,l,k) = elec_coord_tmp(i,l,ipm)
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

  SOFT_TOUCH elec_coord_full_dmc psi_value psi_grad_psi_inv_x psi_grad_psi_inv_y psi_grad_psi_inv_z elec_coord 

  first_loop = .False.

 enddo ! while loop

 double precision :: factor
 factor = 1.d0/block_weight
 SOFT_TOUCH block_weight

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   $X_dmc_block_walk   *= factor
   $X_2_dmc_block_walk *= factor
 endif
"""
for p in properties:
 print (t.replace("$X",p[1]))
END_SHELL

 deallocate(elec_coord_tmp, psi_grad_psi_inv_save, psi_grad_psi_inv_save_tmp)

 do k=1,min(walk_num,walk_num_dmc)
  do l=1,3
   do i=1,elec_num+1
    elec_coord_full(i,l,k) = elec_coord_full_dmc(i,l,k)
   enddo
  enddo
 enddo
 do k=walk_num_dmc+1,walk_num
  do l=1,3
   do i=1,elec_num+1
    elec_coord_full(i,l,k) = elec_coord_full_dmc(i,l,mod(k,walk_num_dmc)+1)
   enddo
  enddo
 enddo
 SOFT_TOUCH elec_coord_full

END_PROVIDER


BEGIN_PROVIDER [ double precision, E_ref ]
  implicit none
  BEGIN_DOC  
!  Weight  of the DMC population
  END_DOC
  E_ref = 0.d0
  call get_simulation_E_ref(E_ref)
END_PROVIDER

 BEGIN_PROVIDER [ integer, trapped_walk, (walk_num_dmc_max) ]
&BEGIN_PROVIDER [ integer,  trapped_walk_max ]
 implicit none
 BEGIN_DOC  
! Number of steps when the walkers were trapped
 END_DOC
 trapped_walk = 0
 trapped_walk_max = 20
END_PROVIDER


BEGIN_PROVIDER [ integer, walk_num_dmc ]
 implicit none
 BEGIN_DOC
 ! Current number of walkers in DMC
 END_DOC
 walk_num_dmc = walk_num
END_PROVIDER

BEGIN_PROVIDER [ integer, walk_num_dmc_max ]
 implicit none
 BEGIN_DOC
 ! Max number of walkers in DMC
 END_DOC
 walk_num_dmc_max = max(3 * walk_num, 30)
END_PROVIDER


BEGIN_PROVIDER [ real, elec_coord_full_dmc, (elec_num+1,3,walk_num_dmc_max)]
 implicit none
 BEGIN_DOC
 ! DMC population
 END_DOC
 integer :: i,k,l
 do k=1,walk_num
   do l=1,3
     do i=1,elec_num+1
       elec_coord_full_dmc(i,l,k) = elec_coord_full(i,l,k)
      enddo
    enddo
  enddo

END_PROVIDER

