subroutine brownian_step(p,q,accepted,delta_x)

  implicit none
  include '../types.F'

  double precision,intent(out)  :: p,q
  logical, intent(out)  :: accepted
  real,intent(out)  :: delta_x

  real              :: xold_x(elec_num+1)
  real              :: xold_y(elec_num+1)
  real              :: xold_z(elec_num+1)
  double precision  :: psiold
  double precision  :: bold0_x(elec_num),bold_x(elec_num),bnew_x(elec_num)
  double precision  :: bold0_y(elec_num),bold_y(elec_num),bnew_y(elec_num)
  double precision  :: bold0_z(elec_num),bold_z(elec_num),bnew_z(elec_num)
  double precision  :: E_old
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xold_x
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xold_y
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xold_z
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: bold_x
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: bold_y
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: bold_z

  integer :: i,l

  psiold = psi_value
  if (psiold == 0.) then
    call abrt(irp_here,'Walker is on the nodal surface.')
  endif

  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(100)
  do i=1,elec_num+1
    xold_x(i) = elec_coord(i,1)
    xold_y(i) = elec_coord(i,2)
    xold_z(i) = elec_coord(i,3)
  enddo
  !DIR$ VECTOR ALIGNED
  bold0_x = psi_grad_psi_inv_x
  bold_x = psi_grad_psi_inv_x
  !DIR$ VECTOR ALIGNED
  bold0_y = psi_grad_psi_inv_y
  bold_y = psi_grad_psi_inv_y
  !DIR$ VECTOR ALIGNED
  bold0_z = psi_grad_psi_inv_z
  bold_z = psi_grad_psi_inv_z

! real :: time_step_inv
! time_step_inv = 1./time_step
! !DIR$ VECTOR ALIGNED
! do i=1,elec_num
!   if (bold_x(i) > time_step_inv) then
!     bold_x(i) = time_step_inv
!   else if (bold_x(i) < -time_step_inv) then
!     bold_x(i) = -time_step_inv
!   endif
!   if (bold_y(i) > time_step_inv) then
!     bold_y(i) = time_step_inv
!   else if (bold_y(i) < -time_step_inv) then
!     bold_y(i) = -time_step_inv
!   endif
!   if (bold_z(i) > time_step_inv) then
!     bold_z(i) = time_step_inv
!   else if (bold_z(i) < -time_step_inv) then
!     bold_z(i) = -time_step_inv
!   endif
! enddo
  

  double precision :: b2old, b2max, tmp
  b2old = 0.d0
  b2max = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(100)
  do i=1,elec_num
    b2old += bold_x(i)*bold_x(i) + bold_y(i)*bold_y(i) + bold_z(i)*bold_z(i)
  enddo

  real :: xdiff_x (elec_num)
  real :: xdiff_y (elec_num)
  real :: xdiff_z (elec_num)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xdiff_x
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xdiff_y
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xdiff_z
  double precision :: gauss 
  do l=1,3
   do i=1,elec_num
    xbrown(i,l) = gauss()*time_step_sq
   enddo
  enddo
  delta_x = 0.
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(100)
  do i=1,elec_num
    xdiff_x(i) = bold_x(i)*time_step + xbrown(i,1)
    xdiff_y(i) = bold_y(i)*time_step + xbrown(i,2)
    xdiff_z(i) = bold_z(i)*time_step + xbrown(i,3)
    delta_x += xdiff_x(i)*xdiff_x(i) + xdiff_y(i)*xdiff_y(i) + xdiff_z(i)*xdiff_z(i)
  enddo
  delta_x = sqrt(delta_x)

  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(100)
  do i=1,elec_num
    elec_coord(i,1) = xold_x(i) + xdiff_x(i)
    elec_coord(i,2) = xold_y(i) + xdiff_y(i)
    elec_coord(i,3) = xold_z(i) + xdiff_z(i)
  enddo

  call update_qmckl_coord()
  TOUCH elec_coord
  
  !DIR$ VECTOR ALIGNED
  bnew_x = psi_grad_psi_inv_x
  !DIR$ VECTOR ALIGNED
  bnew_y = psi_grad_psi_inv_y
  !DIR$ VECTOR ALIGNED
  bnew_z = psi_grad_psi_inv_z

! !DIR$ VECTOR ALIGNED
! do i=1,elec_num
!   if (bnew_x(i) > time_step_inv) then
!     bnew_x(i) = time_step_inv
!   else if (bnew_x(i) < -time_step_inv) then
!     bnew_x(i) = -time_step_inv
!   endif
!   if (bnew_y(i) > time_step_inv) then
!     bnew_y(i) = time_step_inv
!   else if (bnew_y(i) < -time_step_inv) then
!     bnew_y(i) = -time_step_inv
!   endif
!   if (bnew_z(i) > time_step_inv) then
!     bnew_z(i) = time_step_inv
!   else if (bnew_z(i) < -time_step_inv) then
!     bnew_z(i) = -time_step_inv
!   endif
! enddo
  double precision :: ratio
  
  ratio = psi_value/psiold
  ratio *= ratio

  double precision :: b2new
  b2new = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(100)
  do i=1,elec_num
    b2new += bnew_x(i)*bnew_x(i) + bnew_y(i)*bnew_y(i) + bnew_z(i)*bnew_z(i)
  enddo

  double precision :: prod
  prod  = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(100)
  do i=1,elec_num
    prod += ( bnew_x(i)+bold_x(i) )*xdiff_x(i) &
          + ( bnew_y(i)+bold_y(i) )*xdiff_y(i) &
          + ( bnew_z(i)+bold_z(i) )*xdiff_z(i)
  enddo
  
  double precision :: argexpo
  argexpo = 0.5d0*(b2new-b2old)*time_step+prod
  
  p = min (1.d0,ratio*exp(-argexpo))
  q = 1.d0 - p
  
  double precision  :: qmc_ranf
  accepted = p > qmc_ranf()
  if (accepted) then
    accepted_num += 1.
  else
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT(100)
    do i=1,elec_num
     elec_coord(i,1) = xold_x(i)
     elec_coord(i,2) = xold_y(i)
     elec_coord(i,3) = xold_z(i)
    enddo
    !DIR$ VECTOR ALIGNED
    psi_grad_psi_inv_x = bold0_x
    !DIR$ VECTOR ALIGNED
    psi_grad_psi_inv_y = bold0_y
    !DIR$ VECTOR ALIGNED
    psi_grad_psi_inv_z = bold0_z
    psi_value = psiold
    rejected_num += 1.
    call update_qmckl_coord()
    SOFT_TOUCH elec_coord psi_grad_psi_inv_x psi_grad_psi_inv_y psi_grad_psi_inv_z psi_value 
  endif

end
