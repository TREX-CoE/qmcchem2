BEGIN_PROVIDER [ integer, i_state ]
 implicit none
 BEGIN_DOC
 ! Current state
 END_DOC
 i_state = 1
END_PROVIDER

BEGIN_PROVIDER [ integer, N_int ]
 use trexio
 implicit none
 BEGIN_DOC
 ! Number of 64-bit integers needed to represent determinants as binary strings
 END_DOC
 integer(trexio_exit_code) :: rc
 integer :: nint8
 if (use_trexio) then
   rc = trexio_get_int64_num(trexio_file, nint8)
   if (rc /= TREXIO_SUCCESS) stop 'N_int: trexio_get_int64_num failed'
   if (nint8 /= 1) stop 'N_int: int_num /= 1 (mo_num > 64) not supported'
   N_int = nint8
 else
   call get_spindeterminants_n_int(N_int)
 endif
END_PROVIDER

BEGIN_PROVIDER [ integer, bit_kind ]
 implicit none
 BEGIN_DOC
 ! Number of octets per integer storing determinants
 END_DOC
 if (use_trexio) then
   bit_kind = 8
 else
   call get_spindeterminants_bit_kind(bit_kind)
 endif
 ASSERT (bit_kind == 8)
END_PROVIDER

BEGIN_PROVIDER [ integer, N_states ]
 implicit none
 BEGIN_DOC
 ! Number of states in EZFIO file
 END_DOC
 if (use_trexio) then
   N_states = 1
 else
   call get_spindeterminants_n_states(N_states)
 endif
END_PROVIDER

BEGIN_PROVIDER [ integer, det_num_input ]
 use trexio
  implicit none
  BEGIN_DOC
  ! Number of Det_a x Det_b products in input file
  END_DOC
  integer(trexio_exit_code) :: rc
  integer*8 :: ndet8
  if (use_trexio) then
    rc = trexio_read_determinant_num_64(trexio_file, ndet8)
    if (rc /= TREXIO_SUCCESS) stop 'det_num_input: trexio_read_determinant_num_64 failed'
    if (ndet8 < 1_8 .or. ndet8 > int(huge(1))) stop 'det_num_input: bad determinant count'
    det_num_input = int(ndet8)
  else
    call get_spindeterminants_n_det(det_num_input)
  endif
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
  integer*8, allocatable         :: dummy_l(:), dummy_r(:)
  integer :: nalpha, nbeta
  if (use_trexio) then
    allocate (dummy_l(det_num_input), dummy_r(det_num_input))
    call wf_trexio_load(nalpha, nbeta, det_coef_matrix_rows, &
         det_coef_matrix_columns, det_right_coef_matrix_values, &
         dummy_l, dummy_r)
    deallocate (dummy_l, dummy_r)
  else
    allocate (buffer(det_num_input,N_states))
    call get_spindeterminants_psi_coef_matrix_rows(det_coef_matrix_rows)
    call get_spindeterminants_psi_coef_matrix_columns(det_coef_matrix_columns)
    call get_spindeterminants_psi_coef_matrix_values(buffer)
    det_right_coef_matrix_values(:) = buffer(:,i_state)
    deallocate(buffer)
  endif
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
  integer,   allocatable :: dummy_r(:), dummy_c(:)
  integer*8, allocatable :: dummy_ul(:), dummy_ur(:)
  double precision, allocatable :: dummy_v(:)
  if (use_trexio) then
    allocate (dummy_r(det_num_input), dummy_c(det_num_input), &
         dummy_v(det_num_input), dummy_ul(det_num_input), &
         dummy_ur(det_num_input))
    call wf_trexio_load(det_alpha_num, det_beta_num, &
         dummy_r, dummy_c, dummy_v, dummy_ul, dummy_ur)
    deallocate (dummy_r, dummy_c, dummy_v, dummy_ul, dummy_ur)
  else
    call get_spindeterminants_n_det_alpha(det_alpha_num)
    call get_spindeterminants_n_det_beta(det_beta_num)
  endif
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
 integer,   allocatable :: dummy_r(:), dummy_c(:)
 integer*8, allocatable :: dummy_ur(:)
 double precision, allocatable :: dummy_v(:)
 integer :: nalpha, nbeta
 if (use_trexio) then
   if (N_int /= 1) stop 'psi_det_alpha: N_int /= 1 not supported in trexio det path'
   allocate (dummy_r(det_num_input), dummy_c(det_num_input), &
        dummy_v(det_num_input), dummy_ur(det_num_input))
   call wf_trexio_load(nalpha, nbeta, dummy_r, dummy_c, dummy_v, &
        psi_det_alpha(1,:), dummy_ur)
   deallocate (dummy_r, dummy_c, dummy_v, dummy_ur)
 else
   call get_spindeterminants_psi_det_alpha(psi_det_alpha)
 endif
END_PROVIDER

BEGIN_PROVIDER [ integer*8, psi_det_beta, (N_int,det_beta_num) ]
 implicit none
 BEGIN_DOC
 ! Beta determinants
 END_DOC
 integer,   allocatable :: dummy_r(:), dummy_c(:)
 integer*8, allocatable :: dummy_ul(:)
 double precision, allocatable :: dummy_v(:)
 integer :: nalpha, nbeta
 if (use_trexio) then
   if (N_int /= 1) stop 'psi_det_beta: N_int /= 1 not supported in trexio det path'
   allocate (dummy_r(det_num_input), dummy_c(det_num_input), &
        dummy_v(det_num_input), dummy_ul(det_num_input))
   call wf_trexio_load(nalpha, nbeta, dummy_r, dummy_c, dummy_v, &
        dummy_ul, psi_det_beta(1,:))
   deallocate (dummy_r, dummy_c, dummy_v, dummy_ul)
 else
   call get_spindeterminants_psi_det_beta(psi_det_beta)
 endif
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

! ============================================================================
! TREXIO (h5) wavefunction loader — used in TT mode (use_trexio=T).
! Reads the CIPSI determinant wavefunction from the trexio file's
! /determinant group instead of the flat spindeterminants/ files.
! Deterministic: dataserver and every worker compute identical arrays.
! Hard stop on any TREXIO error (no silent fallback to flat).
! The per-process cache (save locals) is built once, on the first call;
! every call then copies the cached arrays into the caller's buffers.
! ============================================================================

! All seven arguments are mandatory: the provider bodies are generated by
! IRPF90 without an explicit interface (it forbids interface blocks in
! provider bodies), so every call passes all seven actuals positionally
! and the arrays are assumed-size dummies (an implicit interface forbids
! allocatable dummies, so the caller pre-allocates every actual to the
! cache size: rows/cols/vals to det_num_input, uniq_l/uniq_r to at
! least nalpha/nbeta — the provider arrays have exactly that size).
! The kinds must match the actuals (no interface to check it): rows/
! cols are integer, vals double precision, uniq_l/uniq_r integer*8.
subroutine wf_trexio_load(nalpha, nbeta, rows, cols, vals, uniq_l, uniq_r)
  use trexio
  implicit none
  integer,            intent(out) :: nalpha, nbeta
  integer,            intent(out) :: rows(*), cols(*)
  double precision,   intent(out) :: vals(*)
  integer*8,          intent(out) :: uniq_l(*), uniq_r(*)
  ! --- per-process cache (built once on the first call) ---
  logical,              save :: wf_tx_loaded = .false.
  integer,              save :: wf_tx_nalpha = 0, wf_tx_nbeta = 0
  integer*8, save, allocatable :: wf_tx_uniq_l(:), wf_tx_uniq_r(:)
  integer,   save, allocatable :: wf_tx_rows(:), wf_tx_cols(:)
  double precision, save, allocatable :: wf_tx_vals(:)
  ! --- locals ---
  integer(trexio_exit_code) :: rc
  integer                   :: int_num
  integer*8                 :: ndet, off, want, got, got2, L, R
  integer*8, allocatable    :: lbuf(:)
  double precision, allocatable :: cbuf(:)
  integer, allocatable      :: tbl_l(:), tbl_r(:)
  integer                   :: tblsize, mask, pos, nuniq_l, nuniq_r
  integer, parameter        :: CH = 262144
  integer*8, parameter      :: CH8 = 262144_8
  integer                   :: k, ktot
  integer*8, external       :: detmix64
  character*(128)           :: msg
  character*(16)            :: s1, s2

  if (.not.wf_tx_loaded) then
    call provide_trexio_file
    if (.not.use_trexio) stop 'wf_trexio_load: called with use_trexio=F'

    rc = trexio_get_int64_num(trexio_file, int_num)
    if (rc /= TREXIO_SUCCESS) stop 'wf_trexio_load: trexio_get_int64_num failed'
    if (int_num /= 1) stop 'wf_trexio_load: int_num /= 1 (mo_num > 64) not supported'

    rc = trexio_read_determinant_num_64(trexio_file, ndet)
    if (rc /= TREXIO_SUCCESS) stop 'wf_trexio_load: trexio_read_determinant_num_64 failed'
    if (ndet < 1_8) stop 'wf_trexio_load: no determinants in trexio file'

    allocate (wf_tx_uniq_l(ndet), wf_tx_uniq_r(ndet))
    allocate (wf_tx_rows(ndet), wf_tx_cols(ndet), wf_tx_vals(ndet))
    nuniq_l = 0
    nuniq_r = 0
    tblsize = 2
    do while (tblsize < 2*int(ndet))
      tblsize = tblsize*2
    enddo
    mask = tblsize - 1
    allocate (tbl_l(tblsize), tbl_r(tblsize))
    tbl_l = 0
    tbl_r = 0
    allocate (lbuf(2*CH), cbuf(CH))

    off = 0_8
    ktot = 0
    do while (off < ndet)
      want = min(CH8, ndet - off)
      got = want
      rc = trexio_read_determinant_list(trexio_file, off, got, lbuf)
      if (rc /= TREXIO_SUCCESS .and. rc /= TREXIO_END) &
        stop 'wf_trexio_load: trexio_read_determinant_list failed'
      got2 = want
      rc = trexio_read_determinant_coefficient(trexio_file, off, got2, cbuf)
      if (rc /= TREXIO_SUCCESS .and. rc /= TREXIO_END) &
        stop 'wf_trexio_load: trexio_read_determinant_coefficient failed'
      if (got < 1_8) stop 'wf_trexio_load: empty determinant chunk'
      if (got2 /= got) stop 'wf_trexio_load: list/coef chunk sizes differ'
      do k = 1, int(got)
        ktot = ktot + 1
        L = lbuf(2*(k-1)+1)   ! int_num == 1: per-det block [left][right]
        R = lbuf(2*(k-1)+2)
        ! --- insert L into left hash table (open addressing, linear probing)
        ! 1-based slot index: hash gives 0..mask, shift to 1..tblsize;
        ! wrap pos=tblsize -> 1 (iand(pos,mask)+1)
        pos = int(iand(detmix64(L), int(mask,8))) + 1
        do while (tbl_l(pos) > 0 .and. wf_tx_uniq_l(tbl_l(pos)) /= L)
          pos = iand(pos, mask) + 1
        enddo
        if (tbl_l(pos) == 0) then
          nuniq_l = nuniq_l + 1
          wf_tx_uniq_l(nuniq_l) = L
          tbl_l(pos) = nuniq_l
        endif
        wf_tx_rows(ktot) = tbl_l(pos)
        ! --- insert R (mirror of the above)
        pos = int(iand(detmix64(R), int(mask,8))) + 1
        do while (tbl_r(pos) > 0 .and. wf_tx_uniq_r(tbl_r(pos)) /= R)
          pos = iand(pos, mask) + 1
        enddo
        if (tbl_r(pos) == 0) then
          nuniq_r = nuniq_r + 1
          wf_tx_uniq_r(nuniq_r) = R
          tbl_r(pos) = nuniq_r
        endif
        wf_tx_cols(ktot) = tbl_r(pos)
        wf_tx_vals(ktot) = cbuf(k)
      enddo
      off = off + got
      if (got < want) exit   ! short chunk (TREXIO_END) -> done
    enddo

    deallocate (tbl_l, tbl_r, lbuf, cbuf)

    if (ktot /= int(ndet)) stop 'wf_trexio_load: determinant count mismatch'
    if (nuniq_l < 1 .or. nuniq_r < 1) stop 'wf_trexio_load: no unique determinants'
    if (abs(wf_tx_vals(1)) < 0.1d0) stop 'wf_trexio_load: leading CIPSI coefficient not dominant'
    wf_tx_nalpha = nuniq_l
    wf_tx_nbeta  = nuniq_r
    wf_tx_loaded = .true.
    write (s1,'(i0)') nuniq_l
    write (s2,'(i0)') nuniq_r
    write (msg,'(a,i0)') 'wf_trexio_load: read ', ndet
    print *, trim(msg) // ' determinants from trexio (nalpha=' // &
          trim(s1) // ', nbeta=' // trim(s2) // ')'
  endif

  nalpha = wf_tx_nalpha
  nbeta  = wf_tx_nbeta
  rows(1:size(wf_tx_rows)) = wf_tx_rows
  cols(1:size(wf_tx_cols)) = wf_tx_cols
  vals(1:size(wf_tx_vals)) = wf_tx_vals
  uniq_l(1:wf_tx_nalpha) = wf_tx_uniq_l(1:wf_tx_nalpha)
  uniq_r(1:wf_tx_nbeta) = wf_tx_uniq_r(1:wf_tx_nbeta)
end

integer*8 function detmix64(x)
! Deterministic 64-bit bit-mix (xorshift triple; no overflow, standards-safe)
  implicit none
  integer*8, intent(in) :: x
  integer*8 :: a
  a = ieor(x, ishft(x,-17_8))
  a = ieor(a, ishft(a,-29_8))
  a = ieor(a, ishft(a,  9_8))
  detmix64 = a
end


