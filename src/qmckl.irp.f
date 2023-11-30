 BEGIN_PROVIDER [ integer*8, qmckl_ctx ]
&BEGIN_PROVIDER [ integer*8, qmckl_mo_num ]

  use qmckl
  implicit none
  BEGIN_DOC
  ! Context for the QMCKL library
  END_DOC
  include 'types.F'

  integer, save :: ifirst = 0
  integer(qmckl_exit_code) :: rc
  integer :: keep(mo_tot_num), i

  if (ifirst == 0) then
    ifirst = 1

    qmckl_ctx = qmckl_context_create()
    
    rc = qmckl_trexio_read(qmckl_ctx, trexio_filename, 1_8*len(trim(trexio_filename)))
    call check_qmckl(rc, irp_here, qmckl_ctx)

    rc = qmckl_set_numprec_precision(qmckl_ctx, precision_bits)

    keep = 0
    do i=1,num_present_mos
      keep(present_mos(i)) = 1
    enddo

    rc = qmckl_mo_basis_select_mo(qmckl_ctx, keep, int(mo_tot_num,8))
    call check_qmckl(rc, irp_here, qmckl_ctx)

    if (do_nucl_fitcusp) then
      rc = qmckl_set_mo_basis_r_cusp(qmckl_ctx, dble(nucl_fitcusp_radius(:)), int(nucl_num,8))
      call check_qmckl(rc, irp_here, qmckl_ctx)
    endif

    rc = qmckl_get_mo_basis_mo_num(qmckl_ctx, qmckl_mo_num)
    call check_qmckl(rc, irp_here, qmckl_ctx)

    if (qmckl_mo_num /= num_present_mos) then
      stop 'qmckl_mo_num /= num_present_mos'
    endif

    ! Jastrow parameters
    if (jast_type == t_Qmckl .or. jpsi_type == t_Qmckl) then

       rc =  qmckl_set_jastrow_champ_spin_independent(qmckl_ctx, 1)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_type_nucl_num(qmckl_ctx, 1_8*jast_qmckl_type_nucl_num)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_type_nucl_vector(qmckl_ctx, 1_8*jast_qmckl_type_nucl_vector-1_8, 1_8*nucl_num)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_rescale_factor_ee(qmckl_ctx, jast_qmckl_rescale_ee)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_rescale_factor_en(qmckl_ctx, jast_qmckl_rescale_en, 1_8*jast_qmckl_type_nucl_num)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_aord_num(qmckl_ctx, jast_qmckl_aord_num*1_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_a_vector(qmckl_ctx, jast_qmckl_a_vector, 1_8*size(jast_qmckl_a_vector))
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_bord_num(qmckl_ctx, jast_qmckl_bord_num*1_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_b_vector(qmckl_ctx, jast_qmckl_b_vector, 1_8*size(jast_qmckl_b_vector))
       call check_qmckl(rc, irp_here, qmckl_ctx)


       rc =  qmckl_set_jastrow_champ_cord_num(qmckl_ctx, jast_qmckl_cord_num*1_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       if (jast_qmckl_cord_num > 0) then
         rc =  qmckl_set_jastrow_champ_c_vector(qmckl_ctx, jast_qmckl_c_vector, 1_8*jast_qmckl_c_vector_size)
         call check_qmckl(rc, irp_here, qmckl_ctx)
       endif

    end if

  end if

END_PROVIDER

subroutine check_qmckl(rc, here, ctx)
  use qmckl
  implicit none
  integer(qmckl_exit_code), intent(inout) :: rc
  character*(*), intent(in) :: here
  integer(qmckl_context), intent(in) :: ctx
  if (rc == QMCKL_SUCCESS) return
  print *, here
  rc = qmckl_check(ctx, rc)
  stop -1
end

BEGIN_PROVIDER [ integer, qmckl_precision ]
 use qmckl
 implicit none
 BEGIN_DOC
 ! Precision used in QMCKL
 END_DOC
 qmckl_precision = qmckl_get_numprec_precision(qmckl_ctx)
END_PROVIDER

BEGIN_PROVIDER [ logical, use_qmckl ]
 implicit none
 BEGIN_DOC
 ! Is true, use TREXIO file
 END_DOC
 use_qmckl = .False.
 call get_simulation_use_qmckl(use_qmckl)
END_PROVIDER

BEGIN_PROVIDER [ double precision, qmckl_ao_vgl, (ao_num, 5, elec_num) ]
 use qmckl
 implicit none
 BEGIN_DOC
 ! AO value, gradients, Laplacian from QMCkl
 END_DOC
 PROVIDE elec_coord

 integer(qmckl_exit_code) :: rc
 rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx, qmckl_ao_vgl, elec_num*ao_num*5_8)
 call check_qmckl(rc, irp_here, qmckl_ctx)

 integer :: i,j
 do j=1,elec_num
   do i=1,5
     call dset_order(qmckl_ao_vgl(1,i,j),ao_nucl_sort_idx, ao_num)
   enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, qmckl_mo_vgl, (qmckl_mo_num, 5, elec_num) ]
 use qmckl
 implicit none
 BEGIN_DOC
 ! MO value, gradients, Laplacian from QMCkl
 END_DOC
 PROVIDE elec_coord
 integer(qmckl_exit_code) :: rc
 rc = qmckl_get_mo_basis_mo_vgl_inplace(qmckl_ctx, qmckl_mo_vgl, walk_num*elec_num*qmckl_mo_num*5_8)
 call check_qmckl(rc, irp_here, qmckl_ctx)

END_PROVIDER


BEGIN_PROVIDER [ double precision, qmckl_mo_value, (qmckl_mo_num, elec_num) ]
 use qmckl
 implicit none
 BEGIN_DOC
 ! MO value, gradients, Laplacian from QMCkl
 END_DOC
 PROVIDE elec_coord
 integer(qmckl_exit_code) :: rc
 rc = qmckl_get_mo_basis_mo_value_inplace(qmckl_ctx, qmckl_mo_value, walk_num*elec_num*qmckl_mo_num)
 call check_qmckl(rc, irp_here, qmckl_ctx)

END_PROVIDER

subroutine update_qmckl_coord()
  implicit none
  double precision :: buffer(elec_num,3)
  integer(qmckl_exit_code) :: rc

  buffer(1:elec_num,1:3) = elec_coord(1:elec_num,1:3)
  rc = qmckl_set_electron_coord(qmckl_ctx, 'T', 1_8, buffer,  3_8*elec_num)
  call check_qmckl(rc, irp_here, qmckl_ctx)
end
