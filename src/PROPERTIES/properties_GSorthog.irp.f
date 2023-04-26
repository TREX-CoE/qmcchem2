

BEGIN_PROVIDER [ double precision, vec_gs_Psi, (size_vec_gs_Psi) ]

  BEGIN_DOC
  !
  ! vec_gs_Psi(i) = < D(i) e^{-J} / Phi >
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = psi_value_inv
  do k = 1, det_num
    i = det_coef_matrix_rows(   k)
    j = det_coef_matrix_columns(k)
    vec_gs_Psi(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_gs_Psi_min = min(vec_gs_Psi_min, minval(vec_gs_Psi))
  vec_gs_Psi_max = max(vec_gs_Psi_max, maxval(vec_gs_Psi))
  SOFT_TOUCH vec_gs_Psi_min vec_gs_Psi_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, vec_gs_Phi, (size_vec_gs_Phi) ]

  BEGIN_DOC
  !
  ! vec_gs_Phi(i) = < D(i) e^{J} / Phi >
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = psidet_right_inv * Jsimple_value
  do k = 1, det_num
    i = det_coef_matrix_rows(   k)
    j = det_coef_matrix_columns(k)
    vec_gs_Phi(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_gs_Phi_min = min(vec_gs_Phi_min, minval(vec_gs_Phi))
  vec_gs_Phi_max = max(vec_gs_Phi_max, maxval(vec_gs_Phi))
  SOFT_TOUCH vec_gs_Phi_min vec_gs_Phi_max
END_PROVIDER
