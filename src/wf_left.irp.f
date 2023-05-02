
! ---

BEGIN_PROVIDER [ double precision, det_left_coef_matrix_values, (det_num_input) ]

  BEGIN_DOC
  ! det_left_coef_matrix in sparse storage (Coordinate format for sparse BLAS)
  END_DOC

  implicit none
  double precision, allocatable :: buffer(:,:)

  print*, ' reading det_left_coef_matrix_values'
  allocate(buffer(det_num_input,N_states))
  call get_spindeterminants_psi_left_coef_matrix_values(buffer)
  det_left_coef_matrix_values(:) = buffer(:,i_state)
  deallocate(buffer)

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, det_left_coef_matrix_dense, (det_alpha_num, det_beta_num) ]

  BEGIN_DOC
  ! Dense version of det_left_coef_matrix
  END_DOC

  implicit none
  integer :: i, j, k

  print*, ' reading det_left_coef_matrix_dense'
  det_left_coef_matrix_dense = 0.d0
  do k = 1, det_num
    i = det_coef_matrix_rows   (k)
    j = det_coef_matrix_columns(k)
    det_left_coef_matrix_dense(i,j) = det_left_coef_matrix_values(k)
  enddo

END_PROVIDER

! ---

