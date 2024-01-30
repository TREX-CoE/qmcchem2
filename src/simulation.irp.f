BEGIN_PROVIDER [ logical, is_worker ]
 implicit none
 BEGIN_DOC
 ! True if the process is a worker process
 END_DOC
 is_worker = .False.
END_PROVIDER


BEGIN_PROVIDER [ integer, walk_num_tot ]
   implicit none
   BEGIN_DOC
   ! Total number of walkers
   END_DOC
   
   walk_num_tot = 1000
   call get_electrons_elec_walk_num_tot(walk_num_tot)
   walk_num_tot = max(walk_num,walk_num_tot)
   call iinfo(irp_here,'walk_num', walk_num_tot)
   
   if (walk_num_tot <= 0) then
     call abrt(irp_here,'Total number of walkers should be > 0')
   endif
END_PROVIDER

 
 BEGIN_PROVIDER [ integer, walk_num ]
&BEGIN_PROVIDER [ integer, walk_num_8 ]
   implicit none
   BEGIN_DOC
   ! Number of walkers
   END_DOC
   
   walk_num = 100
   call get_electrons_elec_walk_num(walk_num)
   call iinfo(irp_here,'walk_num', walk_num)
   
   if (walk_num <= 0) then
     call abrt(irp_here,'Number of walkers should be > 0')
   endif
   integer                        :: mod_align
   walk_num_8 = mod_align(walk_num)
END_PROVIDER
 

BEGIN_PROVIDER [ logical, do_equilibration ]
   implicit none
   BEGIN_DOC
   ! Equilibrate walkers
   END_DOC
   
   do_equilibration = .True.
   if (.not.do_prepare) then
     call get_simulation_equilibration(do_equilibration)
   endif
   
END_PROVIDER

 
BEGIN_PROVIDER [ logical, do_prepare ]
   implicit none
   BEGIN_DOC
   ! If true, prepare new walkers
   END_DOC
   do_prepare = .False.
   
END_PROVIDER

 
BEGIN_PROVIDER [ double precision, block_time ]
   implicit none
   BEGIN_DOC
   ! Wall time requested to realize one block
   END_DOC
   block_time = 30.d0
   integer :: block_time_int
   call get_simulation_block_time(block_time_int)
   
   if (block_time<= 1) then
     call abrt(irp_here,'Block time should be > 1s')
   endif
   double precision, external :: qmc_ranf
   block_time = dble(block_time_int) + qmc_ranf()
   call dinfo(irp_here,'block_time',block_time)
END_PROVIDER

 
BEGIN_PROVIDER [ integer, stop_time ]
   implicit none
   BEGIN_DOC
   ! Termination condition of the run
   END_DOC
   stop_time = 3600*24
   call get_simulation_stop_time(stop_time)
   call iinfo(irp_here,'stop_time',stop_time)
   
   if (stop_time<= 1) then
     call abrt(irp_here,'Stop time should be > 1s')
   endif
END_PROVIDER

 
BEGIN_PROVIDER [ integer, precision_bits ]
   implicit none
   BEGIN_DOC
   ! Termination condition of the run
   END_DOC
   precision_bits = 24
   call get_simulation_precision(precision_bits)
   call iinfo(irp_here,'precision',precision_bits)
   
   if (precision_bits<= 1) then
     call abrt(irp_here,'Precision should be > 1')
   endif
   if (precision_bits> 53) then
     call abrt(irp_here,'Precision should be <= 53')
   endif
END_PROVIDER

 
 BEGIN_PROVIDER [ real, time_step ]
&BEGIN_PROVIDER [ real, time_step_inv ]
&BEGIN_PROVIDER [ double precision, dtime_step ]
   implicit none
   BEGIN_DOC
   ! time_step : The time step of the random walk
   END_DOC
   
   time_step = 0.0
   call get_simulation_time_step(time_step)
   call rinfo(irp_here,'time_step',time_step)
   
   if (time_step <= 0.) then
     call abrt(irp_here,'Time step should be > 0')
   endif
   dtime_step = dble(time_step)
   time_step_inv = 1./time_step
END_PROVIDER
 
 
 BEGIN_PROVIDER [ double precision, time_step_sq ]
&BEGIN_PROVIDER [ double precision, time_step_exp ]
&BEGIN_PROVIDER [ double precision, time_step_exp_sq ]
&BEGIN_PROVIDER [ double precision, time_step_exp_sq_sq ]
   implicit none
   BEGIN_DOC
   !
   ! time_step_sq : sqrt(time_step)
   !
   ! time_step_exp = exp(-time_step)
   !
   ! time_step_exp_sq = sqrt(time_step_exp)
   !
   ! time_step_exp_sq_sq = sqrt(time_step_exp_sq)
   END_DOC
   
   time_step_sq = sqrt(dble(time_step))
   time_step_exp = exp(-dble(time_step))
   time_step_exp_sq = sqrt(time_step_exp)
   time_step_exp_sq_sq = sqrt(time_step_exp_sq)
END_PROVIDER

 
BEGIN_PROVIDER  [ integer, qmc_method ]
   implicit none
   include 'types.F'
   BEGIN_DOC
   ! qmc_method : Calculation method. Can be t_VMC, t_DMC, t_SRMC, t_FKMC
   END_DOC
   character*(32)                 :: method
   method = types(t_VMC)
   call get_simulation_method(method)
   
   if (method == types(t_VMC)) then
     qmc_method = t_VMC
   else if (method == types(t_DMC)) then
     qmc_method = t_DMC
   else if (method == types(t_SRMC)) then
     qmc_method = t_SRMC
   else if (method == types(t_FKMC)) then
     qmc_method = t_FKMC
   else if (method == types(t_PDMC)) then
     qmc_method = t_PDMC
   else
     call abrt(irp_here, 'Method should be ( VMC | DMC | SRMC | FKMC | PDMC )')
   endif
   
   call cinfo(irp_here,'qmc_method',trim(method))
   
END_PROVIDER

 
BEGIN_PROVIDER [ real, events_num ]
   BEGIN_DOC
   ! Number of Monte Carlo events to average
   END_DOC
   events_num = real(walk_num)*real(step_num)
END_PROVIDER

 
BEGIN_PROVIDER [ integer, walk_i ]
   BEGIN_DOC
   ! Current walker
   END_DOC
   walk_i = 1
END_PROVIDER

 
 BEGIN_PROVIDER [ real, accepted_num ]
&BEGIN_PROVIDER [ real, rejected_num ]
   BEGIN_DOC
   ! Number of accepted steps
   ! Number of rejected steps
   END_DOC
   accepted_num= 0.
   rejected_num= 0.
END_PROVIDER
 

BEGIN_PROVIDER [ logical, save_data ]
   implicit none
   BEGIN_DOC
   ! If true, the updated simulation data is saved for restart.
   END_DOC
   save_data = .True.
   call get_simulation_save_data(save_data)
   call linfo(irp_here,'save_data',save_data)
   
END_PROVIDER

 
real function accep_rate()
   if ((accepted_num+rejected_num) > 0.) then
     accep_rate = accepted_num/(accepted_num+rejected_num)
   else
     accep_rate = 0.
   endif
end

 
logical function first_step()
   first_step = (accepted_num+rejected_num < 1.)
end
 

subroutine accep_reset
   FREE accepted_num
   FREE rejected_num
end

 

 
BEGIN_PROVIDER [ integer, print_level ]
   
   BEGIN_DOC
   ! Level of verbosity for standard output printing
   END_DOC
   print_level = 1
   call get_simulation_print_level(print_level)
   
END_PROVIDER

 
BEGIN_PROVIDER [ character*(64), hostname]
 implicit none
 BEGIN_DOC
 ! Name of the current host
 END_DOC
 call HOSTNM(hostname)
END_PROVIDER


 BEGIN_PROVIDER [ real, nucl_fitcusp_factor ]
&BEGIN_PROVIDER [ logical, do_nucl_fitcusp ]
 implicit none
 BEGIN_DOC
 ! The electron-nucleus cusp fitting is done between 0 and r_c,
 ! where r_c is chosen as nucl_fitcusp_factor * (radius_of_1s AO)
 END_DOC
 nucl_fitcusp_factor = 0.
 call get_simulation_nucl_fitcusp_factor(nucl_fitcusp_factor)
 do_nucl_fitcusp = nucl_fitcusp_factor > 0.
 call rinfo(irp_here,'nucl_fitcusp_factor',nucl_fitcusp_factor)

END_PROVIDER


 
BEGIN_PROVIDER [ integer, vmc_algo ]
   implicit none
   include 'types.F'
   BEGIN_DOC
   ! Type of VMC algorithm: Brownian, MTM or Langevin
   END_DOC
   character*(32)                 :: Sampling
   
   Sampling = types(t_Langevin)
   call get_simulation_sampling(Sampling)
   
   vmc_algo = 0
   if (Sampling == types(t_Brownian)) then
     vmc_algo = t_Brownian
   else if (Sampling == types(t_Langevin)) then
     vmc_algo = t_Langevin
     if (qmc_method == t_DMC) then
       stop 'Langevin incompatible with DMC'
     endif
     if (qmc_method == t_SRMC) then
       stop 'Langevin incompatible with SRMC'
     endif
     if (qmc_method == t_FKMC) then
       stop 'Langevin incompatible with FKMC'
     endif
   else if (Sampling == types(t_MTM)) then
     vmc_algo = t_MTM
   else
     call abrt(irp_here,'Sampling should be (Brownian|Langevin)')
   endif
   call cinfo(irp_here,'Sampling',Sampling)
   
   ASSERT (vmc_algo > 0)
END_PROVIDER
 
 
BEGIN_PROVIDER  [ character*(512), ezfio_filename ]
   implicit none
   BEGIN_DOC
   ! Name of the ezfio file.
   ! Defined in init_ezfio_filename
   END_DOC
   integer                        :: command_argument_count
   if (command_argument_count() == 0) then
     ezfio_filename = 'NOT_SET'
     call ezfio_set_file(ezfio_filename)
   else
     call get_command_argument(1,ezfio_filename)
     if (.not.is_worker) then
       call ezfio_set_file(ezfio_filename)
     endif
   endif
END_PROVIDER

BEGIN_PROVIDER [ character*(128), http_server ]
  implicit none
  BEGIN_DOC
  ! Address of the data server
  END_DOC
  integer                        :: command_argument_count
  if (command_argument_count() > 1) then
    call get_command_argument(2,http_server)
  else
    call get_simulation_http_server(http_server)
  endif
END_PROVIDER


subroutine read_do_run(do_run)
 implicit none
 integer, intent(out) :: do_run
 BEGIN_DOC
 ! Read the do_run variable from the ezfio directory
 END_DOC
 include 'types.F'
 do_run = t_Stopping
 call get_simulation_do_run(do_run)
end



BEGIN_PROVIDER [ character*(32), md5_key ]
   implicit none
   BEGIN_DOC
   ! Digest of the input
   END_DOC

   md5_key = ''
   call get_simulation_md5_key(md5_key)

   if (md5_key == '') then
     call abrt(irp_here,'MD5 key of input is absent')
   endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, E_trial ]
 implicit none
 BEGIN_DOC
 ! Energy of the trial wave function
 END_DOC
 call get_simulation_e_trial(E_trial)
END_PROVIDER

