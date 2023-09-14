subroutine draw_init_points
  use qmckl
  implicit none
  BEGIN_DOC
! Place randomly electron around nuclei
  END_DOC
  integer                        :: iwalk
  logical, allocatable           :: do_elec(:)
  integer                        :: acc_num

  real, allocatable              :: xmin(:,:)

  integer                        :: i, j, k, l, kk

  real                           :: norm
  allocate (do_elec(elec_num), xmin(3,elec_num))
  xmin = -huge(1.)
  norm = 0.
  do i=1,elec_alpha_num
    do j=1,ao_num
      norm += mo_coef_transp(i,j)*mo_coef_transp(i,j)
    enddo
  enddo
  norm = sqrt(norm/float(elec_alpha_num))
  call rinfo( irp_here, 'Norm : ', norm )
  mo_coef_transp = mo_coef_transp/norm

  double precision               :: qmc_ranf
  real                           :: mo_max
  do i=1,elec_alpha_num
    l=1
    xmin(1,i) = mo_coef_transp(i,1)*mo_coef_transp(i,1) - 0.001*qmc_ranf()
    do j=2,ao_num
      xmin(2,i) = mo_coef_transp(i,j)*mo_coef_transp(i,j) - 0.001*qmc_ranf()
      if (xmin(2,i) > xmin(1,i)  ) then
        xmin(1,i) = xmin(2,i)
        l = ao_nucl(j)
      endif
    enddo
    xmin(1,i) = nucl_coord(l,1)
    xmin(2,i) = nucl_coord(l,2)
    xmin(3,i) = nucl_coord(l,3)
  enddo

  call iinfo(irp_here, 'Det num = ', det_num )
  do k=1,elec_beta_num
    i = k+elec_alpha_num
    l=1
    xmin(1,i) = mo_coef_transp(k,1)*mo_coef_transp(k,1) - 0.001*qmc_ranf()
    do j=2,ao_num
      xmin(2,i) = mo_coef_transp(k,j)*mo_coef_transp(k,j) - 0.001*qmc_ranf()
      if (xmin(2,i) > xmin(1,i) ) then
        xmin(1,i) = xmin(2,i)
        l = ao_nucl(j)
      endif
    enddo
    xmin(1,i) = nucl_coord(l,1)
    xmin(2,i) = nucl_coord(l,2)
    xmin(3,i) = nucl_coord(l,3)
  enddo

  call rinfo( irp_here, 'time step =', time_step )
  do iwalk=1,walk_num
    print *, 'Generating initial positions for walker', iwalk
    acc_num = 0
    do_elec = .True.
    integer :: iter
    do iter = 1,10000
      if (acc_num >= elec_num) then
        exit
      endif
      double precision               :: gauss
      real                           :: re_compute
      re_compute = 0.
      do while (re_compute < 1.e-6)
        do i=1,elec_num
          if (do_elec(i)) then
            do l=1,3
              elec_coord(i,l) = xmin(l,i) + 1.5*(0.5-qmc_ranf())
            enddo
          endif
        enddo
        call update_qmckl_coord()
        TOUCH elec_coord
        re_compute = minval(nucl_elec_dist(1:nucl_num,1:elec_num))
      enddo

      do i=1,elec_alpha_num
        if (do_elec(i)) then
          if ( mo_value_transp(i,i)**2 >= qmc_ranf()) then
            acc_num += 1
            do_elec(i) = .False.
          endif
        endif
      enddo

      do i=1,elec_beta_num
        if (do_elec(i+elec_alpha_num)) then
          if ( mo_value_transp(i,i+elec_alpha_num)**2 >= qmc_ranf()) then
            acc_num += 1
            do_elec(i+elec_alpha_num) = .False.
          endif
        endif
      enddo

    enddo

    do l=1,3
      do i=1,elec_num+1
        elec_coord_full(i,l,iwalk) = elec_coord(i,l)
      enddo
    enddo
  enddo
  if (.not.is_worker) then
    call ezfio_set_electrons_elec_coord_pool_size(walk_num)
    call ezfio_set_electrons_elec_coord_pool(elec_coord_full)
  endif
  call update_qmckl_coord()
  SOFT_TOUCH elec_coord elec_coord_full
  deallocate (do_elec, xmin)

end


subroutine run_prepare_walkers
  implicit none
  BEGIN_DOC
! Create starting points for walkers
  END_DOC
  include 'types.F'
  integer                        :: istep, iwalk
  integer                        :: i,j, l

  do iwalk=1,walk_num
    do l=1,3
      do i=1,elec_num+1
        elec_coord(i,l) = elec_coord_full(i,l,iwalk)
      enddo
    enddo
    call update_qmckl_coord()
    TOUCH elec_coord

    double precision               :: qmc_ranf, rcond, lambda
    rcond = 100.d0
    lambda = 1.d0
    do while ( (rcond > 3.d0) .or. (rcond < -3.d0) )
      rcond = 0.
      do i=1,elec_alpha_num
        rcond += log(lambda*abs(mo_value_transp(i,i)))
      enddo
      do i=1,elec_beta_num
        rcond += log(lambda*abs(mo_value_transp(i,elec_alpha_num+i)))
      enddo
      if (rcond > 2.d0) then
        lambda = lambda/(1.d0+.1*qmc_ranf())
      endif
      if (rcond< -2.d0) then
        lambda = lambda*(1.d0+.1*qmc_ranf())
      endif
    enddo
    do i=1,ao_num
      !DIR$ VECTOR ALIGNED
      do j=1,mo_num_8
        mo_coef_transp(j,i) *= lambda
      enddo
    enddo
    TOUCH mo_coef_transp
    print *, 'Starting walker ', iwalk

    do istep=1,1000
      if (psidet_right_value == 0.d0) then
        exit
      endif
      prepare_walkers_t = float(istep)/1000.
      TOUCH prepare_walkers_t
      rcond = log(abs(dble(psidet_right_value)))
      real                           :: factor
      rcond = log(abs(dble(psidet_right_value)))
      integer                        :: icount
      icount = 0
      do while ( (rcond > 10.d0) .or. (rcond < -10.d0) )
        icount += 1
        if (icount == 1000) then
          exit
        endif
        if (rcond > 10.d0) then
          factor = 1./(1.+.10) !*qmc_ranf())
        else if (rcond< -10.d0) then
          factor = 1.+.10 !*qmc_ranf()
        endif
        do j=1,ao_num
          !DIR$ VECTOR ALIGNED
          do i=1,mo_num_8
            mo_coef_transp(i,j) *= factor
          enddo
        enddo
        TOUCH mo_coef_transp

        rcond = log(abs(dble(psidet_right_value)))
      enddo
      double precision               :: p,q
      logical                        :: accepted
      real                           :: delta_x
      accepted = .False.
      do icount=1,100
        if (vmc_algo == t_Brownian) then
          call brownian_step(p,q,accepted,delta_x)
        else if (vmc_algo == t_Langevin) then
          call langevin_step(p,q,accepted,delta_x)
        endif
        if (accepted) then
          exit
        endif
      enddo
    enddo
    do l=1,3
      do i=1,elec_num+1
        elec_coord_full(i,l,iwalk) = elec_coord(i,l)
      enddo
    enddo
    TOUCH elec_coord_full
  enddo

end

