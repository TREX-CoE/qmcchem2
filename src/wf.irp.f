BEGIN_PROVIDER [ integer, i_state ]
 implicit none
 BEGIN_DOC
 ! Current state
 END_DOC
 i_state = 1
END_PROVIDER

BEGIN_PROVIDER [ integer, N_int ]
 implicit none
 BEGIN_DOC
 ! Number of 64-bit integers needed to represent determinants as binary strings
 END_DOC
 call get_spindeterminants_n_int(N_int)
END_PROVIDER

BEGIN_PROVIDER [ integer, bit_kind ]
 implicit none
 BEGIN_DOC
 ! Number of octets per integer storing determinants
 END_DOC
 call get_spindeterminants_bit_kind(bit_kind)
 ASSERT (bit_kind == 8)
END_PROVIDER

BEGIN_PROVIDER [ integer, N_states ]
 implicit none
 BEGIN_DOC
 ! Number of states in EZFIO file
 END_DOC
 call get_spindeterminants_n_states(N_states)
END_PROVIDER

BEGIN_PROVIDER [ integer, det_num_input ]
  implicit none
  BEGIN_DOC
  ! Number of Det_a x Det_b products in input file
  END_DOC
  call get_spindeterminants_n_det(det_num_input)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, det_alpha_norm, (det_alpha_num) ]
&BEGIN_PROVIDER [ double precision, det_beta_norm, (det_beta_num) ]
 implicit none
 BEGIN_DOC
 ! Norm of the alpha and beta spin determinants in the wave function:
 !
 ! ||Da||_i \sum_j C_{ij}**2
 END_DOC

 integer :: i,j,k
 double precision :: f

 det_alpha_norm = 0.d0
 det_beta_norm  = 0.d0
 do k=1,det_num
   i = det_coef_matrix_rows(k)
   j = det_coef_matrix_columns(k)
   f = det_right_coef_matrix_values(k) * det_right_coef_matrix_values(k)
   det_alpha_norm(i) += f
   det_beta_norm(j)  += f
 enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, det_right_coef_matrix_values,  (det_num_input) ]
&BEGIN_PROVIDER [ integer,          det_coef_matrix_rows,    (det_num_input) ]
&BEGIN_PROVIDER [ integer,          det_coef_matrix_columns, (det_num_input) ]
  implicit none
  BEGIN_DOC
  ! det_coef_matrix in sparse storage (Coordinate format for sparse BLAS)
  END_DOC
  double precision, allocatable  :: buffer(:,:)
  allocate (buffer(det_num_input,N_states))
  call get_spindeterminants_psi_coef_matrix_rows(det_coef_matrix_rows)
  call get_spindeterminants_psi_coef_matrix_columns(det_coef_matrix_columns)
  call get_spindeterminants_psi_coef_matrix_values(buffer)
  det_right_coef_matrix_values(:) = buffer(:,i_state)
  deallocate(buffer)
END_PROVIDER

BEGIN_PROVIDER [ double precision, det_right_coef_matrix_dense, (det_alpha_num, det_beta_num) ]
 implicit none
 BEGIN_DOC
 ! Dense version of det_coef_matrix
 END_DOC
 integer :: i,j,k
 det_right_coef_matrix_dense = 0.d0
 do k=1,det_num
   i = det_coef_matrix_rows(k)
   j = det_coef_matrix_columns(k)
   det_right_coef_matrix_dense(i,j) = det_right_coef_matrix_values(k)
 enddo
END_PROVIDER


 BEGIN_PROVIDER [ integer, det_num ]
  implicit none
  BEGIN_DOC
  ! Number of Det_a x Det_b products. The determinant basis set is reduced with
  ! the CI threshold
  END_DOC
  integer                        :: i,j,k,l
  double precision               :: f

  double precision :: d_alpha(det_alpha_num), d_beta (det_beta_num)
  integer :: i_alpha(det_alpha_num), i_beta(det_beta_num)
  integer :: iorder(max(det_alpha_num,det_beta_num))
  integer*8, allocatable ::  psi_det_tmp(:,:)
  double precision :: t, norm

  allocate (psi_det_tmp (N_int,max(det_alpha_num,det_beta_num)))

  if (use_svd) then
    t = -1.d0 
  else
    t = ci_threshold
  endif

  ! Compute the norm of the alpha and beta determinants
  d_alpha = 0.d0
  d_beta  = 0.d0
  do k=1,det_num_input
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    f = det_right_coef_matrix_values(k) * det_right_coef_matrix_values(k)
    d_alpha(i) += f
    d_beta (j) += f
  enddo
  t = min(t, maxval(d_alpha))
  t = min(t, maxval(d_beta))

  ! Reorder alpha determinants
  do i=1,det_alpha_num
    iorder(i) = i
    if (d_alpha(i) < t) then
      i_alpha(i) = det_alpha_num+i
    else
      i_alpha(i) = i
    endif
  enddo
  call isort(i_alpha,iorder,det_alpha_num)

  i=det_alpha_num
  do while (i > 0)
    if (i_alpha(i) <= det_alpha_num) then
      det_alpha_num = i
      exit
    else
      i = i-1
    endif
  enddo

  do i=1,det_alpha_num
    psi_det_tmp(:,i) = psi_det_alpha(:,iorder(i))
    i_alpha(iorder(i)) = i
  enddo
  do i=1,det_alpha_num
    psi_det_alpha(:,i) = psi_det_tmp(:,i)
  enddo

  ! Reorder beta determinants
  do i=1,det_beta_num
    iorder(i) = i
    if (d_beta(i) < t) then
      i_beta(i) = det_beta_num+i
    else
      i_beta(i) = i
    endif
  enddo
  call isort(i_beta,iorder,det_beta_num)


  i=det_beta_num
  do while (i > 0)
    if (i_beta(i) <= det_beta_num) then
      det_beta_num = i
      exit
    else
      i = i-1
    endif
  enddo

  do i=1,det_beta_num
    psi_det_tmp(:,i) = psi_det_beta(:,iorder(i))
    i_beta(iorder(i)) = i
  enddo
  do i=1,det_beta_num
    psi_det_beta(:,i) = psi_det_tmp(:,i)
  enddo
  deallocate(psi_det_tmp)


  ! Apply the threshold to the wave function
  l = 1
  norm = 0.d0
  do k=1,det_num_input
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    det_coef_matrix_rows(l)    = i_alpha(i)
    det_coef_matrix_columns(l) = i_beta(j)
    det_right_coef_matrix_values(l)  = det_right_coef_matrix_values(k)
    if ( (d_alpha(i) >= t).and.(d_beta(j) >= t) ) then
      l = l+1
      norm += det_right_coef_matrix_values(k) * det_right_coef_matrix_values(k)
    endif
  enddo
  det_num = l-1
  norm = 1.d0/dsqrt(norm)
  do k=1,det_num
    det_right_coef_matrix_values(k) *= norm
  enddo

  SOFT_TOUCH det_alpha_num det_beta_num det_right_coef_matrix_values det_coef_matrix_rows det_coef_matrix_columns psi_det_beta psi_det_alpha

END_PROVIDER

 BEGIN_PROVIDER [ integer, det_alpha_num ]
&BEGIN_PROVIDER [ integer, det_beta_num  ]
  implicit none
  BEGIN_DOC
  ! Number of alpha and beta determinants
  END_DOC
  call get_spindeterminants_n_det_alpha(det_alpha_num)
  call get_spindeterminants_n_det_beta(det_beta_num)
END_PROVIDER

 BEGIN_PROVIDER [ integer, det_alpha_num_8 ]
&BEGIN_PROVIDER [ integer, det_beta_num_8  ]
  implicit none
  BEGIN_DOC
  ! Number of alpha and beta determinants
  END_DOC
  integer                        :: mod_align
  det_alpha_num_8 = max(4,mod_align(det_alpha_num)) !
  det_beta_num_8  = max(4,mod_align(det_beta_num))  ! Used in 4x unrolling
END_PROVIDER


BEGIN_PROVIDER [ double precision, ci_threshold ]

  implicit none
  BEGIN_DOC
  ! Threshold on absolute value of the CI coefficients of the wave functioE
  END_DOC
  ci_threshold = 0.d0
  call get_simulation_ci_threshold(ci_threshold)
  call dinfo(irp_here,'ci_threshold',ci_threshold)

END_PROVIDER

BEGIN_PROVIDER [ integer*8, psi_det_alpha, (N_int,det_alpha_num) ]
 implicit none
 BEGIN_DOC
 ! Alpha determinants
 END_DOC
 call get_spindeterminants_psi_det_alpha(psi_det_alpha)
END_PROVIDER

BEGIN_PROVIDER [ integer*8, psi_det_beta, (N_int,det_beta_num) ]
 implicit none
 BEGIN_DOC
 ! Beta determinants
 END_DOC
 call get_spindeterminants_psi_det_beta(psi_det_beta)
END_PROVIDER

BEGIN_PROVIDER  [ integer, present_mos, (mo_tot_num) ]
&BEGIN_PROVIDER [ integer, num_present_mos ]
&BEGIN_PROVIDER [ integer, num_present_mos_8 ]
&BEGIN_PROVIDER [ integer, mo_closed_num ]
  implicit none
  BEGIN_DOC
  ! List of used MOs to build the wf in the CI expansion
  END_DOC
  integer*8                      :: tmp_det(N_int)
  integer                        :: i,k
  integer, external              :: mod_align
  PROVIDE det_num

  num_present_mos = mo_tot_num
  do i=1,mo_tot_num
    present_mos(i) = i
  enddo

!---
  present_mos = 0
  tmp_det = 0_8
  do i=1,det_alpha_num
    do k=1,N_int
      tmp_det(k) = ior(tmp_det(k),psi_det_alpha(k,i))
    enddo
  enddo
  do i=1,det_beta_num
    do k=1,N_int
      tmp_det(k) = ior(tmp_det(k),psi_det_beta(k,i))
    enddo
  enddo
  call bitstring_to_list(tmp_det,present_mos,num_present_mos,N_int)
!---

  num_present_mos_8 = mod_align(num_present_mos)

  integer :: list(mo_tot_num), n
  logical :: good

  list = present_mos
  mo_closed_num = elec_beta_num
  do n=1,elec_beta_num
    call list_to_bitstring(tmp_det,present_mos,n,N_int)
    do k=1,N_int
      if (tmp_det(k) == 0_8) then
        exit
      endif
      good = .True.
      do i=1,det_alpha_num
        if (iand(tmp_det(k),psi_det_alpha(k,i)) /= tmp_det(k)) then
          good = .False.
          exit
        endif
      enddo
      if (good) then
        do i=1,det_beta_num
          if (iand(tmp_det(k),psi_det_beta(k,i)) /= tmp_det(k)) then
            good = .False.
            exit
          endif
        enddo
      endif
      if (.not.good) then
        exit
      endif
    enddo
    if (.not.good) then
      mo_closed_num = n-1
      exit
    endif
  enddo
END_PROVIDER

subroutine list_to_bitstring( string, list, n_elements, Nint)
  implicit none
  BEGIN_DOC
  ! Returns the physical string "string(N_int,2)" from the array of
  ! occupations "list(N_int*64,2)
  END_DOC
  integer, intent(in)            :: Nint
  integer*8, intent(out)         :: string(Nint)
  integer, intent(in)            :: list(Nint*64)
  integer, intent(in)            :: n_elements


  integer                        :: i, j
  integer                        :: ipos, iint

  !
  !                                       <== ipos ==>
  !                                                  |
  !                                                  v
  !string :|------------------------|-------------------------|------------------------|
  !        <====      64       ====> <====      64       ====> <====      64       ====>
  !        {        iint            } {         iint         } {         iint         }
  !

  string = 0_8

  do i=1,n_elements
    iint = shiftr(list(i)-1,6) + 1
    ipos = list(i)-shiftl((iint-1),6)-1
    string(iint) = ibset( string(iint), ipos )
  enddo

end


BEGIN_PROVIDER [ integer, det_alpha_order, (det_alpha_num) ]
 implicit none
 BEGIN_DOC
 ! Order in which to compute the alhpa determinants
 END_DOC
 integer :: i
 do i=1,det_alpha_num
   det_alpha_order(i) = i
 enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, det_beta_order, (det_beta_num) ]
 implicit none
 BEGIN_DOC
 ! Order in which to compute the beta determinants
 END_DOC
 integer :: i
 do i=1,det_beta_num
   det_beta_order(i) = i
 enddo
END_PROVIDER

! ---


