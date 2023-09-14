BEGIN_PROVIDER  [ double precision, elec_mass ]
 implicit none
 BEGIN_DOC
! Electron mass in the Langevin algorithm
 END_DOC

 integer :: i
 elec_mass=0.d0
 do i=1,nucl_num
  elec_mass=max(dble(nucl_charge(i)),dble(elec_mass))
 enddo
 elec_mass=elec_mass**1.5
END_PROVIDER

BEGIN_PROVIDER  [ real, elec_impuls_full, (elec_num,4,walk_num)]
 implicit none
 BEGIN_DOC
! Impulsions used in the Langevin algorithm
 END_DOC
 integer :: i,k,l
 double precision :: temp(elec_num,3)

 do k=1,walk_num
  call gauss_array(elec_num*3,temp)
  do l=1,3
   do i=1,elec_num
    elec_impuls_full(i,l,k) = time_step*temp(i,l)
   enddo
  enddo
 enddo
 
END_PROVIDER

BEGIN_PROVIDER  [ real, elec_impuls, (elec_num_8,3)]
 implicit none
 BEGIN_DOC
! Impulsions used in the Langevin algorithm for the current walker
 END_DOC
 integer :: i,l
 do l=1,3
  do i=1,elec_num
   elec_impuls(i,l) = elec_impuls_full(i,l,walk_i)
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER  [ double precision, langevin_sigma, (2) ]
&BEGIN_PROVIDER [ double precision, langevin_c12 ]
&BEGIN_PROVIDER [ double precision, langevin_coef, (4) ]
 implicit none
 BEGIN_DOC
! Quantities needed for the Langevin algorithm.
 END_DOC
 langevin_sigma(1) = sqrt (time_step/elec_mass * &
    (2.d0- (3.d0-4.d0*time_step_exp + time_step_exp**2)/time_step) )
 langevin_sigma(2) = sqrt (elec_mass * (1.d0-time_step_exp**2))
 langevin_c12 = (1.d0-time_step_exp)**2/(langevin_sigma(1)*langevin_sigma(2))
 langevin_coef(1) = 1.d0/elec_mass*time_step*time_step_exp_sq
 langevin_coef(2) = time_step_exp_sq_sq*time_step**2/elec_mass
 langevin_coef(3) = langevin_sigma(2)*sqrt(1.d0-langevin_c12**2)
 langevin_coef(4) = langevin_c12*langevin_sigma(2)/langevin_sigma(1)
END_PROVIDER



subroutine langevin_step(p,q,accepted,delta_x)

  implicit none

  double precision, intent(out)  :: p,q
  logical, intent(out) :: accepted
  real, intent(out) :: delta_x

  real              :: xold(elec_num_8,3)
  real              :: pold(elec_num_8,3)
  double precision  :: psiold
  double precision  :: bold(elec_num_8,3)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xold, bold, pold

  integer :: i,l

  psiold = psi_value
  if (psiold == 0.d0) then
    call abrt(irp_here,'Walker is on the nodal surface.')
  endif

  do l=1,3
   !DIR$ LOOP COUNT (200)
   !DIR$ VECTOR ALIGNED
   do i=1,elec_num
    xold(i,l) = elec_coord(i,l)
    pold(i,l) = elec_impuls(i,l)
   enddo
  enddo
  !DIR$ LOOP COUNT (200)
  !DIR$ VECTOR ALIGNED
  do i=1,elec_num
    bold(i,1) = psi_grad_psi_inv_x(i)
    bold(i,2) = psi_grad_psi_inv_y(i)
    bold(i,3) = psi_grad_psi_inv_z(i)
  enddo

! First move
! W1
  call gauss_array(size(xbrown),xbrown)
  do l=1,3
   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (200)
   do i=1,elec_num
    xbrown(i,l) = xbrown(i,l)*langevin_sigma(1)
   enddo
  enddo

  real :: xdiff (elec_num_8,3)
! q(n+1)
  delta_x = 0.
  do l=1,3
   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (200)
   do i=1,elec_num
    xdiff(i,l) = pold(i,l)*langevin_coef(1) &
               + bold(i,l)*langevin_coef(2) &
               + xbrown(i,l)
    delta_x += xdiff(i,l)*xdiff(i,l)
   enddo
  enddo
  delta_x = sqrt(delta_x)

  do l=1,3
   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (200)
   do i=1,elec_num
    elec_coord(i,l) = elec_coord(i,l) + xdiff(i,l)
   enddo
  enddo

  call update_qmckl_coord()
  TOUCH elec_coord
  
! Second move
! W2
  double precision :: gauss 
  do l=1,3
   !DIR$ VECTOR ALIGNED
   !DIR$ LOOP COUNT (200)
   do i=1,elec_num
      xbrown(i,l) = langevin_coef(3)*gauss() + &
                    langevin_coef(4)*xbrown(i,l)
   enddo
  enddo

! p(n+1)
  !DIR$ LOOP COUNT (200)
  do i=1,elec_num
     elec_impuls(i,1) = time_step_exp*pold(i,1) &
        + time_step*(bold(i,1)+psi_grad_psi_inv_x(i))*time_step_exp_sq &
        + xbrown(i,1)
     elec_impuls(i,2) = time_step_exp*pold(i,2) &
        + time_step*(bold(i,2)+psi_grad_psi_inv_y(i))*time_step_exp_sq &
        + xbrown(i,2)
     elec_impuls(i,3) = time_step_exp*pold(i,3) &
        + time_step*(bold(i,3)+psi_grad_psi_inv_z(i))*time_step_exp_sq &
        + xbrown(i,3)
  enddo
! TOUCH elec_impuls
! (touch moved after acceptation)

  p = 1.d0

  double precision :: ratio
  ratio = (psi_value/psiold)**2

  double precision :: dt_dm
  dt_dm = time_step/elec_mass 
  !DIR$ LOOP COUNT (200)
  do i=1,elec_num

    double precision :: d11
    double precision :: d12
    double precision :: temp
    double precision :: d21
    double precision :: d22
    double precision :: re

    d11 = -xdiff(i,1)  &
          + dt_dm*elec_impuls(i,1)*time_step_exp_sq &
          - time_step*dt_dm*psi_grad_psi_inv_x(i)*time_step_exp_sq_sq

    d12 = xdiff(i,1) &
          - dt_dm*pold(i,1)*time_step_exp_sq  &
          - time_step*dt_dm*bold(i,1)*time_step_exp_sq_sq

    temp = time_step*(psi_grad_psi_inv_x(i)+bold(i,1))*time_step_exp_sq

    d21 = -pold(i,1)+elec_impuls(i,1)*time_step_exp - temp
     

    d22 = elec_impuls(i,1)-pold(i,1)*time_step_exp -temp

    re = (d11/langevin_sigma(1))**2+(d21/langevin_sigma(2))**2  &
     - 2.d0*langevin_c12*d11*d21/(langevin_sigma(1)*langevin_sigma(2))  

    re -= (d12/langevin_sigma(1))**2+(d22/langevin_sigma(2))**2 &
     - 2.d0*langevin_c12*d12*d22/(langevin_sigma(1)*langevin_sigma(2))  

    re = re / (2.d0*(1.d0-langevin_c12**2))

    if (re < 35.d0) then
      p=p*exp(-re)
    else
      p = 0.d0
    endif


    d11 = -xdiff(i,2)  &
          + dt_dm*elec_impuls(i,2)*time_step_exp_sq &
          - time_step*dt_dm*psi_grad_psi_inv_y(i)*time_step_exp_sq_sq

    d12 = xdiff(i,2) &
          - dt_dm*pold(i,2)*time_step_exp_sq  &
          - time_step*dt_dm*bold(i,2)*time_step_exp_sq_sq

    temp = time_step*(psi_grad_psi_inv_y(i)+bold(i,2))*time_step_exp_sq

    d21 = -pold(i,2)+elec_impuls(i,2)*time_step_exp - temp
     

    d22 = elec_impuls(i,2)-pold(i,2)*time_step_exp -temp

    re = (d11/langevin_sigma(1))**2+(d21/langevin_sigma(2))**2  &
     - 2.d0*langevin_c12*d11*d21/(langevin_sigma(1)*langevin_sigma(2))  

    re -= (d12/langevin_sigma(1))**2+(d22/langevin_sigma(2))**2 &
     - 2.d0*langevin_c12*d12*d22/(langevin_sigma(1)*langevin_sigma(2))  

    re = re / (2.d0*(1.d0-langevin_c12**2))

    if (re < 35.d0) then
      p=p*exp(-re)
    else
      p = 0.d0
    endif

    d11 = -xdiff(i,3)  &
          + dt_dm*elec_impuls(i,3)*time_step_exp_sq &
          - time_step*dt_dm*psi_grad_psi_inv_z(i)*time_step_exp_sq_sq

    d12 = xdiff(i,3) &
          - dt_dm*pold(i,3)*time_step_exp_sq  &
          - time_step*dt_dm*bold(i,3)*time_step_exp_sq_sq

    temp = time_step*(psi_grad_psi_inv_z(i)+bold(i,3))*time_step_exp_sq

    d21 = -pold(i,3)+elec_impuls(i,3)*time_step_exp - temp
     

    d22 = elec_impuls(i,3)-pold(i,3)*time_step_exp -temp

    re = (d11/langevin_sigma(1))**2+(d21/langevin_sigma(2))**2  &
     - 2.d0*langevin_c12*d11*d21/(langevin_sigma(1)*langevin_sigma(2))  

    re -= (d12/langevin_sigma(1))**2+(d22/langevin_sigma(2))**2 &
     - 2.d0*langevin_c12*d12*d22/(langevin_sigma(1)*langevin_sigma(2))  

    re = re / (2.d0*(1.d0-langevin_c12**2))

    if (re < 35.d0) then
      p=p*exp(-re)
    else
      p = 0.d0
    endif


   double precision :: sumup, sumdn
   sumup = elec_impuls(i,1)*elec_impuls(i,1) &
          + elec_impuls(i,2)*elec_impuls(i,2) &
          + elec_impuls(i,3)*elec_impuls(i,3)
   sumdn = pold(i,1)*pold(i,1) + pold(i,2)*pold(i,2) + pold(i,3)*pold(i,3)
   re = exp(-(sumup-sumdn)/(2.d0*elec_mass))
   p = p*re

  enddo


  p = min (1.d0,ratio*p)
  q = 1.d0 - p
  
  double precision  :: qmc_ranf
  accepted = p > qmc_ranf()
  if (accepted) then
    accepted_num += 1.
    SOFT_TOUCH elec_impuls
  else
    do l=1,3
     !DIR$ VECTOR ALIGNED
     !DIR$ LOOP COUNT(200)
     do i=1,elec_num
      elec_coord(i,l) = xold(i,l)
      elec_impuls(i,l) = -pold(i,l)
     enddo
    enddo
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT(200)
    do i=1,elec_num
      psi_grad_psi_inv_x(i) = bold(i,1)
      psi_grad_psi_inv_y(i) = bold(i,2)
      psi_grad_psi_inv_z(i) = bold(i,3)
    enddo
    psi_value = psiold
    rejected_num += 1.
    call update_qmckl_coord()
    SOFT_TOUCH elec_coord psi_grad_psi_inv_x psi_grad_psi_inv_y psi_grad_psi_inv_z psi_value 
  endif

end

