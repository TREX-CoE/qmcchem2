program qmcchem_info
  implicit none
  PROVIDE ezfio_filename
  double precision               :: cpu0, cpu1
  character*(8)                  :: str_n
  integer                        :: iargc
  integer                        :: imax
  if (command_argument_count() > 1) then
    call get_command_argument(2,str_n)
    read(str_n,*) imax
  else
    imax = 100
  endif
  print *,  'Using qmckl                              : ', use_qmckl
  print *,  'Precision                                : ', precision_bits
  print *,  'Number of determinants                   : ', det_num
  print *,  'Number of unique alpha/beta determinants : ', det_alpha_num, det_beta_num
  if (use_svd) then
    print *,  'SVD rank                                 : ', n_svd_coefs
  endif
  print *,  'Closed-shell MOs                         : ', mo_closed_num
  print *,  'Number of MOs in determinants            : ', num_present_mos
!  print *,  'Det alpha norm:'
!  print *,  det_alpha_norm
!  print *,  'Det beta norm:'
!  print *,  det_beta_norm
  call step1
  call cpu_time (cpu0)
  call step2(imax)
  call cpu_time (cpu1)
  print *,  'Time for the calculation of E_loc (ms)   : ', 1000.*(cpu1-cpu0)/float(imax)
end

subroutine step1
  implicit none
  print *,  'E_loc                                    : ', E_loc
  PROVIDE E_loc
  if (calc_vec_hjsimple_m_hmu_no3b) then
     print *, 'computing Vec_hjsimple_m_hmu_no3b'
     PROVIDE Vec_hjsimple_m_hmu_no3b 
  end if
end

subroutine step2(imax)
  implicit none
  integer, intent(in)            :: imax
  integer                        :: i
  do i=1,imax
    PROVIDE E_loc
    if (calc_vec_hjsimple_m_hmu_no3b) then
      PROVIDE Vec_hjsimple_m_hmu_no3b 
    end if

    TOUCH elec_coord
  enddo
end
