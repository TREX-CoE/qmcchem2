!==========================================================================!
! PROPERTIES                                                               !
!==========================================================================!


BEGIN_PROVIDER [ double precision, dipole, (3) ]
  implicit none
  BEGIN_DOC
  ! Dipole moment
  END_DOC
  
  integer                        :: i
  dipole = 0.d0
  do i=1,nucl_num
    dipole(1) += nucl_coord(i,1)*nucl_charge(i)
    dipole(2) += nucl_coord(i,2)*nucl_charge(i)
    dipole(3) += nucl_coord(i,3)*nucl_charge(i)
  enddo
  
  do i=1,elec_num
    dipole(1) -= elec_coord(i,1)
    dipole(2) -= elec_coord(i,2)
    dipole(3) -= elec_coord(i,3)
  enddo
  
  dipole *= 2.541765d0
  dipole_min = min(minval(dipole),dipole_min)
  dipole_max = max(minval(dipole),dipole_max)
  SOFT_TOUCH dipole_min dipole_max
END_PROVIDER

BEGIN_PROVIDER  [ double precision, wf_extension ]
  implicit none
  
  BEGIN_DOC
  ! Wave function extension
  END_DOC
  
  wf_extension = 0.d0
  integer                        :: i
  do i=1,elec_num
    wf_extension += elec_coord(i,1)*elec_coord(i,1) + elec_coord(i,2)*elec_coord(i,2) + elec_coord(i,3)*elec_coord(i,3)
  enddo
  
  wf_extension_min = min(wf_extension,wf_extension_min)
  wf_extension_max = max(wf_extension,wf_extension_max)
  SOFT_TOUCH wf_extension_min wf_extension_max
END_PROVIDER

BEGIN_PROVIDER [ double precision, pop_weight ]
 implicit none
 BEGIN_DOC
 ! Weight of the SRMC population
 END_DOC
 include '../types.F'
 if (qmc_method == t_SRMC) then
   pop_weight = srmc_pop_weight_mult
 else if (qmc_method == t_PDMC) then
   pop_weight = pdmc_pop_weight_mult(pdmc_n_diag)
 endif
 pop_weight_min = min(pop_weight,pop_weight_min)
 pop_weight_max = max(pop_weight,pop_weight_max)
 SOFT_TOUCH pop_weight_min pop_weight_max

END_PROVIDER

BEGIN_PROVIDER [ double precision, drift_mod, (size_drift_mod) ]
  implicit none
  BEGIN_DOC
  ! Modulus of the drift per electron
  !
  ! Dimensions : elec_num
  END_DOC
  
  integer                        :: i, j
  do i=1,elec_num
    drift_mod(i) = sqrt(                                             &
        psi_grad_psi_inv_x(i)*psi_grad_psi_inv_x(i) +                &
        psi_grad_psi_inv_y(i)*psi_grad_psi_inv_y(i) +                &
        psi_grad_psi_inv_z(i)*psi_grad_psi_inv_z(i) )
  enddo
  drift_mod_min = min(minval(drift_mod),drift_mod_min)
  drift_mod_max = max(maxval(drift_mod),drift_mod_max)
  SOFT_TOUCH drift_mod_min drift_mod_max
  
END_PROVIDER


