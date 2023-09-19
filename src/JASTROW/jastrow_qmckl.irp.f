! QMCkl Jastrow
! --------------

BEGIN_PROVIDER [ integer, jast_qmckl_type_nucl_num ]
 implicit none
 BEGIN_DOC
 ! Number of different nuclei types in QMCkl jastrow
 END_DOC

 jast_qmckl_type_nucl_num = -1
 call get_jastrow_jast_qmckl_type_nucl_num(jast_qmckl_type_nucl_num)
 if (jast_qmckl_type_nucl_num <= 0) then
    call abrt(irp_here, 'jast_qmckl_type_nucl_num <= 0')
 endif

END_PROVIDER

BEGIN_PROVIDER [ integer, jast_qmckl_type_nucl_vector, (nucl_num) ]
 implicit none
 BEGIN_DOC
 ! Nucleus type in QMCkl jastrow
 END_DOC

 jast_qmckl_type_nucl_vector = -1
 call get_jastrow_jast_qmckl_type_nucl_vector(jast_qmckl_type_nucl_vector)

 integer                        :: i
 do i=1,nucl_num
   if (jast_qmckl_type_nucl_vector(i) <= 0) then
     call abrt(irp_here,'jast_qmckl_type_nucl_vector(i) <= 0')
   endif
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, jast_qmckl_rescale_ee ]
 implicit none
 BEGIN_DOC
 ! Rescaling factor for electron-electron in QMCkl Jastrow
 END_DOC

 jast_qmckl_rescale_ee = 0.d0
 call get_jastrow_jast_qmckl_rescale_ee(jast_qmckl_rescale_ee)
 if (jast_qmckl_rescale_ee <= 0.d0) then
    call abrt(irp_here, 'jast_qmckl_rescale_ee <= 0.d0')
 endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, jast_qmckl_rescale_en, (jast_qmckl_type_nucl_num) ]
 implicit none
 BEGIN_DOC
 ! Rescaling factor for electron-nucleus in QMCkl Jastrow
 END_DOC

 jast_qmckl_rescale_en = 0.d0
 call get_jastrow_jast_qmckl_rescale_en(jast_qmckl_rescale_en)

 integer                        :: i
 do i=1,jast_qmckl_type_nucl_num
   if (jast_qmckl_rescale_en(i) <= 0.d0) then
     call abrt(irp_here,'jast_qmckl_rescale_en(i) <= 0.d0')
   endif
 enddo

END_PROVIDER

BEGIN_PROVIDER [ integer, jast_qmckl_aord_num ]
 implicit none
 BEGIN_DOC
 ! Order of polynomials in e-n parameters of QMCkl jastrow
 END_DOC

 jast_qmckl_aord_num = 0
 call get_jastrow_jast_qmckl_aord_num(jast_qmckl_aord_num)
 if (jast_qmckl_aord_num < 0) then
    call abrt(irp_here, 'jast_qmckl_aord_num < 0')
 endif

END_PROVIDER

BEGIN_PROVIDER [ integer, jast_qmckl_bord_num ]
 implicit none
 BEGIN_DOC
 ! Order of polynomials in e-e parameters of QMCkl jastrow
 END_DOC

 jast_qmckl_bord_num = 0
 call get_jastrow_jast_qmckl_bord_num(jast_qmckl_bord_num)
 if (jast_qmckl_bord_num < 0) then
    call abrt(irp_here, 'jast_qmckl_bord_num < 0')
 endif

END_PROVIDER

BEGIN_PROVIDER [ integer, jast_qmckl_cord_num ]
 implicit none
 BEGIN_DOC
 ! Order of polynomials in e-e-n parameters of QMCkl jastrow
 END_DOC

 jast_qmckl_cord_num = 0
 call get_jastrow_jast_qmckl_cord_num(jast_qmckl_cord_num)
 if (jast_qmckl_cord_num < 0) then
    call abrt(irp_here, 'jast_qmckl_cord_num < 0')
 endif

END_PROVIDER



BEGIN_PROVIDER [ integer, jast_qmckl_c_vector_size ]
 implicit none
 BEGIN_DOC
 ! Number of electron-electron-nucleus parameters in QMCkl Jastrow
 END_DOC

 jast_qmckl_c_vector_size = 1
 call get_jastrow_jast_qmckl_c_vector_size(jast_qmckl_c_vector_size)

END_PROVIDER


BEGIN_PROVIDER [ double precision, jast_qmckl_a_vector, (jast_qmckl_type_nucl_num*(jast_qmckl_aord_num+1))]
 implicit none
 BEGIN_DOC
 ! electron-nucleus parameters in QMCkl Jastrow
 END_DOC

 jast_qmckl_a_vector = 0.d0
 call get_jastrow_jast_qmckl_a_vector(jast_qmckl_a_vector)

END_PROVIDER


BEGIN_PROVIDER [ double precision, jast_qmckl_b_vector, (jast_qmckl_bord_num+1)]
 implicit none
 BEGIN_DOC
 ! electron-electron parameters in QMCkl Jastrow
 END_DOC

 jast_qmckl_b_vector = 0.d0
 call get_jastrow_jast_qmckl_b_vector(jast_qmckl_b_vector)

END_PROVIDER


BEGIN_PROVIDER [ double precision, jast_qmckl_c_vector, (jast_qmckl_c_vector_size)]
 implicit none
 BEGIN_DOC
 ! electron-electron-nucleus parameters in QMCkl Jastrow
 END_DOC

 jast_qmckl_c_vector = 0.d0
 call get_jastrow_jast_qmckl_c_vector(jast_qmckl_c_vector)

END_PROVIDER

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


