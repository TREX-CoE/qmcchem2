program main_prepare_walkers
 implicit none

 ! Delete old walkers files if they exist
 call system ('bash -c "rm -f '//trim(ezfio_filename)//'/electrons/elec_coord_pool{.gz,_size}"')

 call set_parameters
 call draw_init_points
 call run_prepare_walkers
 call save_elec_coord_full
 call run
end

subroutine run
 implicit none
 TOUCH elec_coord_full

 E_loc_min = huge(1.d0)
 E_loc_max = -huge(1.d0)
 print *,  '<E_loc>           min              max' 
 integer :: i
 do i=1,4
   call test_block
   call save_elec_coord_full
   TOUCH elec_coord_full
   E_loc_min = huge(1.d0)
   E_loc_max = -huge(1.d0)
 enddo
 print *,  'Generated ', walk_num, ' walkers'
end

subroutine set_parameters
 implicit none
 include '../types.F'
 do_prepare = .True.
 call ezfio_set_properties_e_loc(.True.)
 qmc_method = t_VMC
 vmc_algo = t_Langevin
 ci_threshold = 0.99999
 time_step = 0.1
 block_time = 2.
 prepare_walkers_t = 0.
 time_step_inv = 1./time_step
 dtime_step = dble(time_step)
 SOFT_TOUCH vmc_algo ci_threshold time_step prepare_walkers_t block_time time_step_inv dtime_step do_prepare
end

subroutine test_block
  implicit none
  print *,  real(E_loc_vmc_block_walk), real(E_loc_min), real(E_loc_max)
end
