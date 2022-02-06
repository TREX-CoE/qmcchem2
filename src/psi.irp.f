 BEGIN_PROVIDER [ double precision, psi_value ]
  implicit none
  BEGIN_DOC
  ! Value of the wave function
  END_DOC
  
  psi_value = psidet_value*jast_value
  if (psi_value == 0.d0) then
    call abrt(irp_here,"Value of the wave function is 0.")
  endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_value_inv ]
  implicit none
  BEGIN_DOC
! 1./psi_value
  END_DOC
  psi_value_inv = 1.d0/psi_value
END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_value_inv2 ]
 implicit none
 BEGIN_DOC
 ! 1./(psi_value)**2
 END_DOC
 psi_value_inv2 = psi_value_inv*psi_value_inv
END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_lapl, (elec_num_8) ]
  implicit none
  BEGIN_DOC
  ! Laplacian of the wave function
  END_DOC
  
  integer                        :: i, j
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (100)
  do j=1,elec_num
    psi_lapl(j) = jast_value*(psidet_grad_lapl(4,j) + psidet_value*jast_lapl_jast_inv(j) + 2.d0*(&
        psidet_grad_lapl(1,j)*jast_grad_jast_inv_x(j) +                   &
        psidet_grad_lapl(2,j)*jast_grad_jast_inv_y(j) +                   &
        psidet_grad_lapl(3,j)*jast_grad_jast_inv_z(j) ))
  enddo
  
END_PROVIDER

 BEGIN_PROVIDER [ double precision, psi_grad_psi_inv_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, psi_grad_psi_inv_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision, psi_grad_psi_inv_z, (elec_num_8) ]
  implicit none
  BEGIN_DOC
! grad(psi)/psi
  END_DOC
  
  integer                        :: j
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (100)
  do j=1,elec_num
    psi_grad_psi_inv_x(j) = psidet_grad_lapl(1,j)*psidet_inv + jast_grad_jast_inv_x(j)
    psi_grad_psi_inv_y(j) = psidet_grad_lapl(2,j)*psidet_inv + jast_grad_jast_inv_y(j)
    psi_grad_psi_inv_z(j) = psidet_grad_lapl(3,j)*psidet_inv + jast_grad_jast_inv_z(j)
  enddo
  
END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_lapl_psi_inv, (elec_num_8) ]
  implicit none
  BEGIN_DOC
  ! (Laplacian psi) / psi
  END_DOC
  
  integer                        :: i, j
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (100)
  do j=1,elec_num
    psi_lapl_psi_inv(j) = psi_lapl(j)*psi_value_inv
  enddo
  
END_PROVIDER

