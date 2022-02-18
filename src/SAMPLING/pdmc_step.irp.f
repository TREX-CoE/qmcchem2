! Providers of *_pdmc_block_walk
!==============================
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

t = """
 BEGIN_PROVIDER [ $T, $X_pdmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_pdmc_block_walk_kahan $D2 ]
&BEGIN_PROVIDER [ $T, $X_2_pdmc_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_2_pdmc_block_walk_kahan $D2 ]
 implicit none
 BEGIN_DOC
 ! PDMC averages of $X. Computed in E_loc_pdmc_block_walk
 END_DOC
 $X_pdmc_block_walk = 0.d0
 $X_pdmc_block_walk_kahan = 0.d0
 $X_2_pdmc_block_walk = 0.d0
 $X_2_pdmc_block_walk_kahan = 0.d0
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



 BEGIN_PROVIDER [ double precision, E_loc_pdmc_block_walk              ]
&BEGIN_PROVIDER [ double precision, E_loc_2_pdmc_block_walk            ]
&BEGIN_PROVIDER [ double precision, E_loc_pdmc_block_walk_kahan  , (3) ]
&BEGIN_PROVIDER [ double precision, E_loc_2_pdmc_block_walk_kahan, (3) ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Properties averaged over the block using the PDMC method
 END_DOC

  real, allocatable :: elec_coord_tmp(:,:,:)
  integer :: mod_align
  double precision :: E_loc_save(4,walk_num_dmc_max)
  double precision :: psi_value_save(walk_num)
  double precision :: psi_value_save_tmp(walk_num)
  double precision :: pdmc_weight(walk_num)
  double precision, allocatable :: psi_grad_psi_inv_save(:,:,:)
  double precision, allocatable :: psi_grad_psi_inv_save_tmp(:,:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_grad_psi_inv_save_tmp
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: E_loc_save
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: psi_value_save_tmp
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: pdmc_weight
  allocate ( psi_grad_psi_inv_save(elec_num_8,3,walk_num) ,          &
       psi_grad_psi_inv_save_tmp(elec_num_8,3,walk_num) ,            &
      elec_coord_tmp(mod_align(elec_num+1),3,walk_num) )
  psi_value_save = 0.d0
  psi_value_save_tmp = 0.d0
  pdmc_weight = 1.d0

! Initialization
 if (vmc_algo /= t_Brownian) then
   call abrt(irp_here,'PDMC should run with Brownian algorithm')
 endif

 integer :: k, i_walk, i_step

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   !DIR$ VECTOR ALIGNED
   $X_pdmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_pdmc_block_walk_kahan = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_pdmc_block_walk = 0.d0
   !DIR$ VECTOR ALIGNED
   $X_2_pdmc_block_walk_kahan = 0.d0
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
 double precision               :: delta, thr

 thr = 2.d0/time_step_sq

 logical :: first_loop
 first_loop = .True.
 if (walk_num > 1) then
   call abrt(irp_here,'walk_num > 1')
 endif

  integer :: info
!  double precision :: H(0:pdmc_n_diag/2,0:pdmc_n_diag/2), S(0:pdmc_n_diag/2,0:pdmc_n_diag/2), w(0:pdmc_n_diag/2), work(3*pdmc_n_diag+1)
!  H = 0.d0
!  S = 0.d0

 do while (loop)

  i_walk = 1

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
      E_loc_save(:,i_walk) = E_loc
    endif

   double precision               :: p,q
   real                           :: delta_x
   logical                        :: accepted
   call brownian_step(p,q,accepted,delta_x)

!   if ( psi_value * psi_value_save(i_walk) >= 0.d0 ) then

!2    delta = (E_loc+E_loc_save(1,i_walk))*0.5d0
!3    delta = (5.d0 * E_loc + 8.d0 * E_loc_save(1,i_walk) - E_loc_save(2,i_walk))/12.d0
      delta = (9.d0*E_loc+19.d0*E_loc_save(1,i_walk)-5.d0*E_loc_save(2,i_walk)+E_loc_save(3,i_walk))/24.d0
!      delta = -((-251.d0*E_loc)-646.d0*E_loc_save(1,i_walk)+264.d0*E_loc_save(2,i_walk)-&
!         106.d0*E_loc_save(3,i_walk)+19.d0*E_loc_save(4,i_walk))/720.d0

     delta = (delta - E_ref)*p

     if (delta >= 0.d0) then
       pdmc_weight(i_walk) = dexp(-dtime_step*delta)
     else
       pdmc_weight(i_walk) = 2.d0-dexp(dtime_step*delta)
     endif
     elec_coord(elec_num+1,1) += p*time_step
     elec_coord(elec_num+1,2)  = E_loc
     elec_coord(elec_num+1,3)  = pdmc_weight(i_walk) * pdmc_pop_weight_mult(pdmc_n_diag)
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
     E_loc_save(4,i_walk) = E_loc_save(3,i_walk)
     E_loc_save(3,i_walk) = E_loc_save(2,i_walk)
     E_loc_save(2,i_walk) = E_loc_save(1,i_walk)
     E_loc_save(1,i_walk) = E_loc


  if (do_print_dmc_data) then
    do k=1,walk_num
      double precision, external :: qmc_ranf
      if (qmc_ranf() < 0.001) then
       print *, '--'
       do i=1,elec_num
         print *, elec_coord_full(i,1:3,k)
       enddo
       print *, 'w=', pdmc_pop_weight_mult(pdmc_n_diag) * pdmc_weight(i_walk)
      endif
    enddo
  endif

     if (dabs(pdmc_weight(i_walk)*pdmc_pop_weight_mult(pdmc_n_diag)) > 1.d-15) then
       dmc_zv_weight = 1.d0/(pdmc_weight(i_walk)*pdmc_pop_weight_mult(pdmc_n_diag))
       dmc_zv_weight_half = 1.d0/(pdmc_weight(i_walk)*pdmc_pop_weight_mult(pdmc_n_diag/2))
     else
       dmc_zv_weight = 0.d0
       dmc_zv_weight_half = 0.d0
     endif
     TOUCH dmc_zv_weight dmc_zv_weight_half

!    do i=1,pdmc_n_diag+1
!      E_loc_zv(i) =   E_loc * pdmc_pop_weight_mult(i-1) * pdmc_weight(i_walk) * dmc_zv_weight + (E_trial-E_loc) * dmc_zv_weight
!      E_loc_zv(i+pdmc_n_diag+1) = pdmc_pop_weight_mult(i-1) * pdmc_weight(i_walk) * dmc_zv_weight
!    enddo

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
     if (calc_$X) then
   ! Kahan's summation algorithm to compute these sums reducing the rounding error:
   !  $X_pdmc_block_walk    += $X * pdmc_pop_weight_mult(pdmc_n_diag) * pdmc_weight(i_walk)
   !  $X_2_pdmc_block_walk  += $X_2 * pdmc_pop_weight_mult(pdmc_n_diag) * pdmc_weight(i_walk)
   ! see http://en.wikipedia.org/wiki/Kahan_summation_algorithm

      $X_pdmc_block_walk_kahan($D2 3) = $X * pdmc_pop_weight_mult(pdmc_n_diag) * pdmc_weight(i_walk) - $X_pdmc_block_walk_kahan($D2 1)
      $X_pdmc_block_walk_kahan($D2 2) = $X_pdmc_block_walk $D1  + $X_pdmc_block_walk_kahan($D2 3)
      $X_pdmc_block_walk_kahan($D2 1) = ($X_pdmc_block_walk_kahan($D2 2) - $X_pdmc_block_walk $D1 ) &
          - $X_pdmc_block_walk_kahan($D2 3)
      $X_pdmc_block_walk $D1  =  $X_pdmc_block_walk_kahan($D2 2)


      $X_2_pdmc_block_walk_kahan($D2 3) = $X_2 * pdmc_pop_weight_mult(pdmc_n_diag) * pdmc_weight(i_walk) - $X_2_pdmc_block_walk_kahan($D2 1)
      $X_2_pdmc_block_walk_kahan($D2 2) = $X_2_pdmc_block_walk $D1 + $X_2_pdmc_block_walk_kahan($D2 3)
      $X_2_pdmc_block_walk_kahan($D2 1) = ($X_2_pdmc_block_walk_kahan($D2 2) - $X_2_pdmc_block_walk $D1 ) &
          - $X_2_pdmc_block_walk_kahan($D2 3)
      $X_2_pdmc_block_walk $D1 =  $X_2_pdmc_block_walk_kahan($D2 2)
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

    block_weight += pdmc_pop_weight_mult(pdmc_n_diag) * pdmc_weight(i_walk)

    pdmc_pop_weight_mult(0) = 1.d0/pdmc_weight(i_walk)
!    do k=0,pdmc_n_diag/2
!      do l=0,pdmc_n_diag/2
!        H(k,l) += E_loc*pdmc_pop_weight_mult(k+l) * pdmc_weight(i_walk)
!        S(k,l) += pdmc_pop_weight_mult(k+l) * pdmc_weight(i_walk)
!      enddo
!    enddo
!    H = H + (E_trial - E_loc)

!  else
!     pdmc_weight(i_walk) = 1.d0
!     pdmc_pop_weight(:,:) = 1.d0
!     pdmc_pop_weight_mult(:) = 1.d0
!  endif

  do k=1,pdmc_n_diag
    ! Move to the next projection step
    if (pdmc_projection(pdmc_n_diag) > 0) then
      pdmc_projection_step(k) = mod(pdmc_projection_step(k),pdmc_projection(k))+1
    else
      pdmc_projection_step(k) = 1
    endif

    ! Eventually, recompute the weight of the population
    if (pdmc_projection_step(k) == k) then
      pdmc_pop_weight_mult(k) = 1.d0
      do l=1,pdmc_projection(k)
        pdmc_pop_weight_mult(k) *= pdmc_pop_weight(l,k)
      enddo
    endif
    ! Remove contribution of the old value of the weight at the new
    ! projection step

    pdmc_pop_weight_mult(k) *= 1.d0/pdmc_pop_weight(pdmc_projection_step(k),k)
    pdmc_pop_weight(pdmc_projection_step(k),k) = pdmc_weight(i_walk)/dble(walk_num)

    ! Update the running population weight
    pdmc_pop_weight_mult(k) *= pdmc_pop_weight(pdmc_projection_step(k),k)

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

  SOFT_TOUCH elec_coord_full pdmc_pop_weight_mult

  first_loop = .False.

 enddo


 double precision :: factor
 factor = 1.d0/block_weight
 SOFT_TOUCH block_weight

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *
t = """
 if (calc_$X) then
   $X_pdmc_block_walk   *= factor
   $X_2_pdmc_block_walk *= factor
 endif
"""
for p in properties:
 print (t.replace("$X",p[1]))
END_SHELL

!  H(0,0) = H(3,3)
!  H(1,0) = H(4,3)
!  H(0,1) = H(3,4)
!  H(1,1) = H(4,4)
!  S(0,0) = S(3,3)
!  S(1,0) = S(4,3)
!  S(0,1) = S(3,4)
!  S(1,1) = S(4,4)
!
!  print *,  H(0,0)/S(0,0)
!  print *,  H(1,1)/S(1,1)
!  print *,  ''
!
!  call dsygv(1, 'N', 'U', pdmc_n_diag/2+1, H, pdmc_n_diag/2+1, S, pdmc_n_diag/2+1, w, work, 3*(pdmc_n_diag+1), info)
!  call dsygv(1, 'N', 'U', 2, H, pdmc_n_diag/2+1, S, pdmc_n_diag/2+1, w, work, 3*(pdmc_n_diag+1), info)
!  E_loc_zv_diag_pdmc_block_walk = w(0)
!  print *,  w

 deallocate ( elec_coord_tmp, psi_grad_psi_inv_save, psi_grad_psi_inv_save_tmp )

END_PROVIDER


 BEGIN_PROVIDER [ integer, pdmc_projection, (pdmc_n_diag) ]
&BEGIN_PROVIDER [ integer, pdmc_projection_step, (pdmc_n_diag) ]
 implicit none
 BEGIN_DOC
! Number of projection steps for PDMC
 END_DOC
 real :: pdmc_projection_time
 pdmc_projection_time = 1.
 call get_simulation_srmc_projection_time(pdmc_projection_time)
 pdmc_projection(pdmc_n_diag) = int( pdmc_projection_time/time_step)
 integer :: k
 do k=1,pdmc_n_diag-1
   pdmc_projection(k) = k*pdmc_projection(pdmc_n_diag)/pdmc_n_diag
 enddo
 pdmc_projection_step(:) = 0
END_PROVIDER

BEGIN_PROVIDER [ double precision, pdmc_pop_weight, (0:pdmc_projection(pdmc_n_diag)+1,pdmc_n_diag) ]
 implicit none
 BEGIN_DOC
! Population weight of PDMC
 END_DOC
 pdmc_pop_weight(:,:) = 1.d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, pdmc_pop_weight_mult, (0:pdmc_n_diag) ]
 implicit none
 BEGIN_DOC
! Population weight of PDMC
 END_DOC
 pdmc_pop_weight_mult(:) = 1.d0
END_PROVIDER


BEGIN_PROVIDER [ integer, pdmc_n_diag ]
 implicit none
 BEGIN_DOC
! Size of the matrix to diagonalize
 END_DOC
 pdmc_n_diag = 8
END_PROVIDER



