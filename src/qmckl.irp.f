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
    qmckl_ctx = qmckl_context_create()
    ifirst = 1
    rc = qmckl_trexio_read(qmckl_ctx, trexio_filename, 1_8*len(trim(trexio_filename)))
    call check_qmckl(rc, irp_here, qmckl_ctx)

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
    if (jast_type == t_Qmckl) then
       rc =  qmckl_set_jastrow_champ_type_nucl_num     (qmckl_ctx, 2_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_type_nucl_vector  (qmckl_ctx, (/0_8,1_8,1_8/), 1_8*nucl_num)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_rescale_factor_ee (qmckl_ctx, 0.6d0)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_rescale_factor_en (qmckl_ctx, (/0.6d0, 0.6d0 /), 2_8 )
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_aord_num          (qmckl_ctx, 5_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_bord_num          (qmckl_ctx, 5_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_cord_num          (qmckl_ctx, 0_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

!      double precision :: a_vector(12) = dble(&
!               (/ 0.00000000,  0.00000000, -0.71168405, -0.44415699, -0.13865109,  0.07002267 , &
!                  0.00000000,  0.00000000, -0.11379992,  0.04542846,  0.01696997, -0.01809299 /) )
!    
!      double precision :: b_vector(6) = dble(&
!               (/  0.00000000,  0.65603311,  0.14581988,  0.03138163,  0.00153156, -0.00447302 /) )
!    
!      double precision :: c_vector(46) = &
!               (/ 1.06384279d0, -1.44303973d0, -0.92409833d0,  0.11845356d0, -0.02980776d0, &
!                  1.07048863d0,  0.06009623d0, -0.01854872d0, -0.00915398d0,  0.01324198d0, &
!                 -0.00504959d0, -0.01202497d0, -0.00531644d0,  0.15101629d0, -0.00723831d0, &
!                 -0.00384182d0, -0.00295036d0, -0.00114583d0,  0.00158107d0, -0.00078107d0, &
!                 -0.00080000d0, -0.14140576d0, -0.00237271d0, -0.03006706d0,  0.01537009d0, &
!                 -0.02327226d0,  0.16502789d0, -0.01458259d0, -0.09946065d0,  0.00850029d0, &
!                 -0.02969361d0, -0.01159547d0,  0.00516313d0,  0.00405247d0, -0.02200886d0, &
!                  0.03376709d0,  0.01277767d0, -0.01523013d0, -0.00739224d0, -0.00463953d0, &
!                  0.00003174d0, -0.01421128d0,  0.00808140d0,  0.00612988d0, -0.00610632d0, &
!                  0.01926215d0 /)

       double precision :: a_vector(12) = dble(&
                (/ 0.00000000 , 0.00000000, -0.45105821, -0.23519218, -0.03825391,  0.10072866, &
                   0.00000000 , 0.00000000, -0.06930592, -0.02909224, -0.00134650,  0.01477242 /) )
     
       double precision :: b_vector(6) = dble(&
                (/  0.00000000,  0.00000000,  0.29217862, -0.00450671, -0.02925982, -0.01381532 /) )
     
       double precision :: c_vector(46)
       c_vector = 0.d0


       rc =  qmckl_set_jastrow_champ_a_vector(qmckl_ctx, a_vector, 12_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

       rc =  qmckl_set_jastrow_champ_b_vector(qmckl_ctx, b_vector, 6_8)
       call check_qmckl(rc, irp_here, qmckl_ctx)

!       rc =  qmckl_set_jastrow_champ_c_vector(qmckl_ctx, c_vector, 46_8)
!       call check_qmckl(rc, irp_here, qmckl_ctx)
    end if

  end if

END_PROVIDER

subroutine check_qmckl(rc, here, ctx)
  use qmckl
  implicit none
  integer(qmckl_exit_code), intent(inout) :: rc
  character*(*), intent(in) :: here
  integer(qmckl_context), intent(in) :: ctx
  character*(128) :: msg
  if (rc == QMCKL_SUCCESS) return
  print *, here
  rc = qmckl_check(ctx, rc)
  print *, msg
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
