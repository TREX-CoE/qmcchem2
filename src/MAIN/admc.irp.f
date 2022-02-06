program admc
  call run
  call ezfio_finish
end program admc

subroutine run
  implicit none

    include '../types.F'

    integer          :: iter
    integer          :: k_full  ! Index of walkers in elec_coord_full
    integer          :: k_new   ! Index of walkers in elec_coord_new
    integer          :: iw      ! Number of copies in branching
    integer          :: l

    real, allocatable :: elec_coord_new(:,:,:)

    double precision :: w
    double precision, allocatable :: E_out(:), w_sum(:)

    double precision, external :: qmc_ranf

    allocate(elec_coord_new(elec_num+1,3,walk_num))
    allocate(E_out(walk_num), w_sum(walk_num))

  ! Initialization
  if (vmc_algo /= t_Brownian) then
     call abrt(irp_here,'DMC should run with Brownian algorithm')
  endif

  do iter=1,1000
     call read_coords()
     k_new  = 1
     do k_full=1,walk_num
        call pdmc_trajectory(k_full, w, E_out(k_full), w_sum(k_full))

          ! Find number of copies
          iw = int(w)
          w = w - int(w)
          if (qmc_ranf() < w) then
             iw = iw+1
          end if

          ! Duplicate walker
          do l=1,iw
            elec_coord_new(1:elec_num+1,1:3,k_new) = &
                 elec_coord(1:elec_num+1,1:3)
            k_new = k_new+1
            if (k_new >= walk_num) exit
          end do


       if (k_new >= walk_num) then
         w_sum(k_full+1:) = 0.d0
         exit
       end if
     end do

     k_new = k_new-1
     elec_coord_full(1:elec_num+1,1:3,1:k_new) = &
          elec_coord_new(1:elec_num+1,1:3,1:k_new)

     call write_coords(k_new)
     call write_energy(walk_num, E_out, w_sum)
  end do

end subroutine run

subroutine read_coords()
  implicit none

  integer :: i, k, l

  do k=1,walk_num
     do l=1,3
       do i=1,elec_num+1
        read(*,*) elec_coord_full(i,l,k)
       end do
     end do
  end do

  SOFT_TOUCH elec_coord_full

end subroutine read_coords
subroutine write_coords(nw)
  implicit none
  integer, intent(in) :: nw

  integer :: i, k, l

  write(*,*) nw
  do k=1,nw
     do l=1,3
       do i=1,elec_num+1
        write(*,*) elec_coord_full(i,l,k)
       end do
     end do
  end do

end subroutine write_coords
subroutine write_energy(walk_num_, E_out, w_sum)
  implicit none
  integer, intent(in) :: walk_num_
  double precision, intent(in) :: E_out(walk_num_)
  double precision, intent(in) :: w_sum(walk_num_)

  integer :: i, k
  double precision :: E, S

  E = 0.d0
  S = 0.d0
  do k=1,walk_num
     S = S + w_sum(k)
     E = E + w_sum(k) * E_out(k)
  end do
  write(*,*) 'E', E/S, S

end subroutine write_energy
subroutine pdmc_trajectory(k_full, pdmc_weight, E_out, w_sum)
  implicit none

  integer, intent(in) :: k_full
  double precision, intent(out) :: pdmc_weight, E_out, w_sum

    integer          :: i,j,l
    double precision :: delta

    ! If true, continue to make more steps
    logical :: loop

    ! Max number of steps
    integer :: imax
    integer, parameter :: nmax=10000

    ! Brownian step variables
    double precision               :: p,q
    real                           :: delta_x
    logical                        :: accepted

    ! Local energies from the past
    double precision :: E_loc_save(4)
    double precision :: w


  elec_coord(1:elec_num+1,1:3) = elec_coord_full(1:elec_num+1,1:3,k_full)
  TOUCH elec_coord

  E_out = 0.d0
  w_sum = 0.d0

  E_loc_save(1:4) = E_loc

  pdmc_weight = 1.d0
  loop = .True.

  do imax = 1, nmax

     call brownian_step(p,q,accepted,delta_x)

!    delta = (9.d0*E_loc+19.d0*E_loc_save(1)-5.d0*E_loc_save(2)+E_loc_save(3))/24.d0
     delta = E_loc
     delta = (delta - E_ref)*p

     if (delta >= 0.d0) then
       w = dexp(-dtime_step*delta)
     else
       w = 2.d0-dexp(dtime_step*delta)
     endif
     pdmc_weight = pdmc_weight * w
     elec_coord(elec_num+1,1) += p*time_step
     elec_coord(elec_num+1,2)  = E_loc
     elec_coord(elec_num+1,3)  = pdmc_weight
     if (accepted) then
        E_loc_save(4) = E_loc_save(3)
        E_loc_save(3) = E_loc_save(2)
        E_loc_save(2) = E_loc_save(1)
        E_loc_save(1) = E_loc
     endif

     w_sum = w_sum + pdmc_weight
     E_out = E_out + pdmc_weight * E_loc

     loop = pdmc_weight > 0.5d0 .and. pdmc_weight < 2.0d0
     if (.not.loop) exit

  end do

  E_out = E_out / w_sum

end subroutine pdmc_trajectory
