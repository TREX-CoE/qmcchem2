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

  jast_elec_Qmckl_value(1:elec_num) = qmckl_jast_value / dble(elec_num)

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
  jast_elec_Qmckl_grad_x(i) = qmckl_jast_grad_lapl(i,1)
  jast_elec_Qmckl_grad_y(i) = qmckl_jast_grad_lapl(i,3)
  jast_elec_Qmckl_grad_z(i) = qmckl_jast_grad_lapl(i,3)
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_lapl, (elec_num_8) ]
 implicit none
 BEGIN_DOC  
! Laplacian of the Jastrow factor
 END_DOC

 integer :: i

 do i=1,elec_num
  jast_elec_Qmckl_lapl(i) = qmckl_jast_grad_lapl(i,4)
 enddo

END_PROVIDER

