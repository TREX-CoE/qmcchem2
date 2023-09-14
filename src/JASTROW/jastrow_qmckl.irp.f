! QMCkl Jastrow
! --------------

! ---

BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_value, (elec_num_8)  ]
 use qmckl
 implicit none
 BEGIN_DOC
! Value of QMCkl Jastrow
 END_DOC

 integer(qmckl_exit_code) :: rc
 integer :: i
 double precision :: v(1)
 PROVIDE elec_coord

 rc = qmckl_get_jastrow_champ_factor_ee(qmckl_ctx, v, 1_8)
 call check_qmckl(rc, irp_here, qmckl_ctx)

 v(1) = v(1)/dble(elec_num)
 do i=1,elec_num
   jast_elec_Qmckl_value(i) = v(1)
 enddo

 rc = qmckl_get_jastrow_champ_factor_en(qmckl_ctx, v, 1_8)
 call check_qmckl(rc, irp_here, qmckl_ctx)

 v(1) = v(1)/dble(elec_num)
 do i=1,elec_num
   jast_elec_Qmckl_value(i) += v(1)
 enddo

 rc = qmckl_get_jastrow_champ_factor_een(qmckl_ctx, v, 1_8)
 call check_qmckl(rc, irp_here, qmckl_ctx)

 v(1) = v(1)/dble(elec_num)
 do i=1,elec_num
   jast_elec_Qmckl_value(i) += v(1)
 enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_grad_x, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_grad_y, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_grad_z, (elec_num_8) ]
&BEGIN_PROVIDER [ double precision , jast_elec_Qmckl_lapl  , (elec_num_8) ]
 use qmckl
 implicit none
 BEGIN_DOC
! Gradient and Laplacian of the QMCkl Jastrow factor
 END_DOC

 integer(qmckl_exit_code) :: rc
 double precision, allocatable :: buffer(:,:)
 integer :: i
 PROVIDE elec_coord

 allocate(buffer(elec_num,4))

 rc = qmckl_get_jastrow_champ_factor_ee_gl(qmckl_ctx, buffer, 4_8*elec_num)
 call check_qmckl(rc, irp_here, qmckl_ctx)

 do i=1,elec_num
   jast_elec_Qmckl_grad_x(i) = buffer(i,1)
   jast_elec_Qmckl_grad_y(i) = buffer(i,2)
   jast_elec_Qmckl_grad_z(i) = buffer(i,3)
   jast_elec_Qmckl_lapl  (i) = buffer(i,4)
 enddo

 rc = qmckl_get_jastrow_champ_factor_en_gl(qmckl_ctx, buffer, 4_8*elec_num)
 call check_qmckl(rc, irp_here, qmckl_ctx)

 do i=1,elec_num
   jast_elec_Qmckl_grad_x(i) += buffer(i,1)
   jast_elec_Qmckl_grad_y(i) += buffer(i,2)
   jast_elec_Qmckl_grad_z(i) += buffer(i,3)
   jast_elec_Qmckl_lapl  (i) += buffer(i,4)
 enddo

 rc = qmckl_get_jastrow_champ_factor_een_gl(qmckl_ctx, buffer, 4_8*elec_num)
 call check_qmckl(rc, irp_here, qmckl_ctx)

 do i=1,elec_num
   jast_elec_Qmckl_grad_x(i) += buffer(i,1)
   jast_elec_Qmckl_grad_y(i) += buffer(i,2)
   jast_elec_Qmckl_grad_z(i) += buffer(i,3)
   jast_elec_Qmckl_lapl  (i) += buffer(i,4)
 enddo
END_PROVIDER


