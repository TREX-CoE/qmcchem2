!==========================================================================!
! DIMENSIONS
!==========================================================================!


BEGIN_PROVIDER [ double precision, single_det_E_kin ]
  implicit none
  BEGIN_DOC
  ! Electronic Kinetic energy : -1/2 (Lapl.Psi)/Psi
  END_DOC
  integer                        :: i
  single_det_E_kin = 0.d0
  do i=1,elec_num
    single_det_E_kin -= 0.5d0*single_det_lapl(i)/single_det_value
  enddo

END_PROVIDER


BEGIN_PROVIDER  [ double precision, single_det_E_loc ]
  implicit none
  BEGIN_DOC
  ! Local energy : single_det_E_kin + E_pot + E_nucl
  END_DOC

  single_det_E_loc = single_det_E_kin + E_pot + E_nucl
END_PROVIDER


BEGIN_PROVIDER [ double precision, E_pot_grad, (elec_num,3) ]
  implicit none
  BEGIN_DOC
  ! Gradient of the Electronic Potential energy
  END_DOC

  integer                        :: i,j
  double precision               :: dinv
  do i=1,elec_num
    E_pot_grad(i,1) = 0.d0
    E_pot_grad(i,2) = 0.d0
    E_pot_grad(i,3) = 0.d0
  enddo
  do j=1,elec_num
    do i=1,j-1
      dinv = elec_dist_inv(i,j)
      dinv = dinv*dinv*dinv
      E_pot_grad(i,1) -= elec_dist_vec_x(i,j)*dinv
      E_pot_grad(i,2) -= elec_dist_vec_y(i,j)*dinv
      E_pot_grad(i,3) -= elec_dist_vec_z(i,j)*dinv
    enddo
    do i=j+1,elec_num
      dinv = elec_dist_inv(i,j)
      dinv = dinv*dinv*dinv
      E_pot_grad(i,1) -= elec_dist_vec_x(i,j)*dinv
      E_pot_grad(i,2) -= elec_dist_vec_y(i,j)*dinv
      E_pot_grad(i,3) -= elec_dist_vec_z(i,j)*dinv
    enddo
  enddo
  do i=1,elec_num
    do j=1,nucl_num
      dinv = nucl_charge(j)*nucl_elec_dist_inv(j,i)**3
      E_pot_grad(i,1) += nucl_elec_dist_vec(1,j,i)*dinv
      E_pot_grad(i,2) += nucl_elec_dist_vec(2,j,i)*dinv
      E_pot_grad(i,3) += nucl_elec_dist_vec(3,j,i)*dinv
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, E_pot_elec, (elec_num) ]
  implicit none
  BEGIN_DOC
  ! Electronic Potential energy
  END_DOC

  integer                        :: i, j
  if (do_pseudo) then
    do i=1,elec_num
      E_pot_elec(i) = v_pseudo_local(i) + pseudo_non_local(i)
    enddo
  else
    do i=1,elec_num
      E_pot_elec(i) = 0.d0
    enddo
  endif

  do i=1,elec_num
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT(50)
    do j=1,elec_num
      E_pot_elec(i) = E_pot_elec(i) + 0.5d0*elec_dist_inv(j,i)
    enddo
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT(50)
    do j=1,nucl_num
      E_pot_elec(i) = E_pot_elec(i) - nucl_charge(j)*nucl_elec_dist_inv(j,i)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, E_pot_elec_one, (elec_num) ]
  implicit none
  BEGIN_DOC
  ! Electronic Potential energy
  END_DOC

  integer                        :: i, j
  do i=1,elec_num
    E_pot_elec_one(i) = 0.d0
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT(100)
    do j=1,nucl_num
      E_pot_elec_one(i) -= nucl_charge(j)*nucl_elec_dist_inv(j,i)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, E_pot_elec_two, (elec_num) ]
  implicit none
  BEGIN_DOC
  ! Electronic Potential energy
  END_DOC

  integer                        :: i, j
  do i=1,elec_num
    E_pot_elec_two(i) = 0.d0
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT(200)
    do j=1,elec_num
      if (j==i) then
        cycle
      endif
      E_pot_elec_two(i) += 0.5d0*elec_dist_inv(j,i)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, E_kin_elec, (elec_num) ]
  implicit none

  BEGIN_DOC
  ! Electronic Kinetic energy : -1/2 (Lapl.Psi)/Psi
  END_DOC

  integer                        :: i
  do i=1,elec_num
    E_kin_elec(i) = -0.5d0*psi_lapl_psi_inv(i)
  enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, dmc_zv_weight ]
 implicit none
 BEGIN_DOC
 ! Weight for Zero-variance in DMC
 END_DOC
 dmc_zv_weight = 1.d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, dmc_zv_weight_half ]
 implicit none
 BEGIN_DOC
 ! Weight for Zero-variance in DMC
 END_DOC
 dmc_zv_weight_half = 1.d0
END_PROVIDER


!==========================================================================!
! PROPERTIES                                                              !
!==========================================================================!

BEGIN_PROVIDER [ double precision, E_nucl ]
  implicit none
  BEGIN_DOC
  ! Nuclear potential energy
  END_DOC

  E_nucl = 0.d0
  integer                        :: i, j
  do i=1,nucl_num
    do j=1,i-1
      E_nucl += nucl_charge(i)*nucl_charge(j)/nucl_dist(j,i)
    enddo
  enddo

  E_nucl_min = min(E_nucl,E_nucl_min)
  E_nucl_max = max(E_nucl,E_nucl_max)
  SOFT_TOUCH E_nucl_min E_nucl_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, E_pot ]
  implicit none
  BEGIN_DOC
  ! Electronic Potential energy
  END_DOC

  E_pot = 0.d0
  integer                        :: i, j
  do i=1,elec_num
    E_pot += E_pot_elec(i)
  enddo

  E_pot_min = min(E_pot,E_pot_min)
  E_pot_max = max(E_pot,E_pot_max)
  SOFT_TOUCH E_pot_min E_pot_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, E_kin ]
  implicit none
  BEGIN_DOC
  ! Electronic Kinetic energy : -1/2 (Lapl.Psi)/Psi
  END_DOC

  E_kin = 0.d0

  integer                        :: i
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(200)
  do i=1,elec_num
    E_kin -= 0.5d0*psi_lapl_psi_inv(i)
  enddo

  E_kin_min = min(E_kin,E_kin_min)
  E_kin_max = max(E_kin,E_kin_max)
  SOFT_TOUCH E_kin_min E_kin_max
END_PROVIDER


BEGIN_PROVIDER  [ double precision, E_loc ]
  implicit none
  include '../types.F'
  BEGIN_DOC
  ! Local energy : E_kin + E_pot + E_nucl
  END_DOC

  integer                        :: i
  E_loc = E_nucl
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT(200)
  do i=1,elec_num
    E_loc += E_kin_elec(i) + E_pot_elec(i)
  enddo

!  ! Avoid divergence of E_loc and population explosion
!  if (do_pseudo .and. (.not. (qmc_method == t_VMC) )) then
!    double precision :: delta_e
!    delta_e = E_loc-E_ref
!    E_loc = E_ref + delta_e * dexp(-dabs(delta_e)*time_step)
!  endif
  E_loc_min = min(E_loc,E_loc_min)
  E_loc_max = max(E_loc,E_loc_max)
  SOFT_TOUCH E_loc_min E_loc_max
END_PROVIDER


!BEGIN_PROVIDER [ double precision, E_loc_zv, ((pdmc_n_diag+1)*2) ]
BEGIN_PROVIDER [ double precision, E_loc_zv ]
 implicit none
 BEGIN_DOC
 ! Zero-variance parameter on E_loc
 END_DOC
   E_loc_zv = E_loc
   E_loc_zv += (E_trial-E_loc) * dmc_zv_weight
!   E_loc_zv += - time_step*(E_trial**2 + 1.44341217940434 - E_loc**2)*dmc_zv_weight
!  E_loc_zv(3)  = dmc_zv_weight_half
!  E_loc_zv(:)  = 0.d0

END_PROVIDER
