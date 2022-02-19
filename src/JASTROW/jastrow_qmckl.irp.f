! Qmckl Jastrow
! --------------

BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_value, (elec_num_8)  ]
  use qmckl
  implicit none
  BEGIN_DOC  
  ! QMCkl Jastrow
  END_DOC
  integer :: i

  PROVIDE qmckl_ctx
  integer :: rc

  double precision :: tmp

!  rc = qmckl_get_jastrow_factor_ee(qmckl_ctx, tmp)
!  if (rc /= QMCKL_SUCCESS) stop -1
!
!  tmp = dlog(tmp)/dble(elec_num)
!  jast_elec_Qmckl_value(1:elec_num) = tmp

END_PROVIDER

 BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_grad_z, (elec_num_8) ]
 implicit none
 BEGIN_DOC  
! Gradient of the Jastrow factor
 END_DOC

 integer :: i

 do i=1,elec_num
  jast_elec_Qmckl_grad_x(i) = 0.d0
  jast_elec_Qmckl_grad_y(i) = 0.d0
  jast_elec_Qmckl_grad_z(i) = 0.d0
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_lapl, (elec_num_8) ]
 implicit none
 BEGIN_DOC  
! Laplacian of the Jastrow factor
 END_DOC

 integer :: i

 do i=1,elec_num
  jast_elec_Qmckl_lapl(i) = 0.d0
 enddo

END_PROVIDER

