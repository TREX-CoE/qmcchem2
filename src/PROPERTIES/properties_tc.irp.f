! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJpsi, (size_vec_HJpsi) ]

  BEGIN_DOC
  !
  ! < D_I Eloc_Jpsi e^{-2 Jpsi} / Phi >
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w
  
  ! dont touch it: jast sampling and not jast_psi
  w = Eloc_Jpsi * jast_value_inv * Psi_value_inv

  do k = 1, det_num
    i = det_coef_matrix_rows   (k)
    j = det_coef_matrix_columns(k)
    vec_HJpsi(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJpsi_min = min(vec_HJpsi_min, minval(vec_HJpsi))
  vec_HJpsi_max = max(vec_HJpsi_max, maxval(vec_HJpsi))
  SOFT_TOUCH vec_HJpsi_min vec_HJpsi_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJsimple, (size_vec_HJsimple) ]

  BEGIN_DOC
  !
  ! < D_I Eloc_Jsimple e^{-2 Jsimple} / Phi >
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w
  
  w = Eloc_Jsimple * jast_value_inv * Psi_value_inv
  do k = 1, det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    vec_HJsimple(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJsimple_min = min(vec_HJsimple_min, minval(vec_HJsimple))
  vec_HJsimple_max = max(vec_HJsimple_max, maxval(vec_HJsimple))
  SOFT_TOUCH vec_HJsimple_min vec_HJsimple_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJpsi_m_H, (size_vec_HJpsi_m_H) ]

  BEGIN_DOC
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  ! dont touch it: jast sampling and not jast_psi
  w = deltaE_Jpsi * Psi_value_inv * jast_value_inv
  do k = 1, det_num
     i = det_coef_matrix_rows   (k)
     j = det_coef_matrix_columns(k)
     vec_HJpsi_m_H(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJpsi_m_H_min = min(vec_HJpsi_m_H_min, minval(vec_HJpsi_m_H))
  vec_HJpsi_m_H_max = max(vec_HJpsi_m_H_max, maxval(vec_HJpsi_m_H))
  SOFT_TOUCH vec_HJpsi_m_H_min vec_HJpsi_m_H_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJmu_m_H, (size_vec_HJmu_m_H) ]

  BEGIN_DOC
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = deltaE_Jmu1b * Psi_value_inv * jast_value_inv
  do k = 1, det_num
     i = det_coef_matrix_rows   (k)
     j = det_coef_matrix_columns(k)
     vec_HJmu_m_H(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJmu_m_H_min = min(vec_HJmu_m_H_min, minval(vec_HJmu_m_H))
  vec_HJmu_m_H_max = max(vec_HJmu_m_H_max, maxval(vec_HJmu_m_H))
  SOFT_TOUCH vec_HJmu_m_H_min vec_HJmu_m_H_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJsimple_m_H, (size_vec_HJsimple_m_H) ]

  BEGIN_DOC
  !
  ! < D_I * deltaE_Jsimple * e^{-2 Jsimple} / Phi > 
  !
  ! with:
  !        deltaE_Jsimple = Eloc_Jsimple - Eloc_noJ
  ! and:
  !        Eloc_Jsimple = H (Phi e^{J}) / (Phi e^{J})
  !        Eloc_noJ     = H (Phi      ) / (Phi      )
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = deltaE_Jsimple * Psi_value_inv * jast_value_inv
  do k = 1, det_num
     i = det_coef_matrix_rows(   k)
     j = det_coef_matrix_columns(k)
     vec_HJsimple_m_H(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJsimple_m_H_min = min(vec_HJsimple_m_H_min, minval(vec_HJsimple_m_H))
  vec_HJsimple_m_H_max = max(vec_HJsimple_m_H_max, maxval(vec_HJsimple_m_H))
  SOFT_TOUCH vec_HJsimple_m_H_min vec_HJsimple_m_H_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJsimple_m_Herf, (size_vec_HJsimple_m_Herf) ]

  BEGIN_DOC
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = EJsimple_m_Eerf * Psi_value_inv * jast_value_inv
  do k = 1, det_num
     i = det_coef_matrix_rows   (k)
     j = det_coef_matrix_columns(k)
     vec_HJsimple_m_Herf(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJsimple_m_Herf_min = min(vec_HJsimple_m_Herf_min, minval(vec_HJsimple_m_Herf))
  vec_HJsimple_m_Herf_max = max(vec_HJsimple_m_Herf_max, maxval(vec_HJsimple_m_Herf))
  SOFT_TOUCH vec_HJsimple_m_Herf_min vec_HJsimple_m_Herf_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJsimple_m_Hmu, (size_vec_HJsimple_m_Hmu) ]

  BEGIN_DOC
  !
  ! < D_I [ Eloc_HJsimple - Eloc_Jmu1b ] * e^{-2J} / Phi >
  ! 
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = EJsimple_m_EJmu1b * Psi_value_inv * jast_value_inv
  do k = 1, det_num
     i = det_coef_matrix_rows(   k)
     j = det_coef_matrix_columns(k)
     vec_HJsimple_m_Hmu(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJsimple_m_Hmu_min = min(vec_HJsimple_m_Hmu_min, minval(vec_HJsimple_m_Hmu))
  vec_HJsimple_m_Hmu_max = max(vec_HJsimple_m_Hmu_max, maxval(vec_HJsimple_m_Hmu))
  SOFT_TOUCH vec_HJsimple_m_Hmu_min vec_HJsimple_m_Hmu_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, vec_HJsimple_m_Hmu_no3b, (size_vec_HJsimple_m_Hmu_no3b) ]

  BEGIN_DOC
  !
  !        Eloc_Jsimple    = H (Phi e^{Jsimple}) / (Phi e^{Jsimple})
  !        Eloc_Jmu1b_no3b = H (Phi e^{Jmu1b}  ) / (Phi e^{Jmu1b}  ) - deltaE_Jmu_grad_3b
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = EJsimple_m_EJmu1b_no3b * Psi_value_inv * jast_value_inv
  do k = 1, det_num
     i = det_coef_matrix_rows(   k)
     j = det_coef_matrix_columns(k)
     vec_HJsimple_m_Hmu_no3b(k) = w * det_right_alpha_value(i) * det_right_beta_value(j) 
  enddo

  vec_HJsimple_m_Hmu_no3b_min = min(vec_HJsimple_m_Hmu_no3b_min, minval(vec_HJsimple_m_Hmu_no3b))
  vec_HJsimple_m_Hmu_no3b_max = max(vec_HJsimple_m_Hmu_no3b_max, maxval(vec_HJsimple_m_Hmu_no3b))
  SOFT_TOUCH vec_HJsimple_m_Hmu_no3b_min vec_HJsimple_m_Hmu_no3b_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, psi_norm ]

  BEGIN_DOC
  ! <1/J^2>
  END_DOC

  implicit none

  psi_norm = jast_value_inv * jast_value_inv

  psi_norm_min = min(psi_norm_min, psi_norm)
  psi_norm_max = max(psi_norm_max, psi_norm)
  SOFT_TOUCH psi_norm_min psi_norm_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

!BEGIN_PROVIDER [ double precision, chi_sq ]
!
!  implicit none
!  integer          :: i, j, k
!  double precision :: w, tmp
!  
!  tmp = EJsimple_m_EJmu1b_no3b * Psi_value_inv * jast_value_inv
!  w   = tmp * tmp
! 
!  chi_sq = 0.d0
!  do k = 1, det_num
!    i = det_coef_matrix_rows(   k)
!    j = det_coef_matrix_columns(k)
!    tmp = det_right_alpha_value(i) * det_right_beta_value(j) 
!    chi_sq += w * tmp * tmp * det_right_coef_matrix_values(k) 
!  enddo
!
!  chi_sq_min = min(chi_sq_min, chi_sq)
!  chi_sq_max = max(chi_sq_max, chi_sq)
!  SOFT_TOUCH chi_sq_min chi_sq_max
!END_PROVIDER

!BEGIN_PROVIDER [ double precision, chi_sq ]
!
!  implicit none
!  
!  chi_sq = EJsimple_m_EJmu1b_no3b * EJsimple_m_EJmu1b_no3b
!
!  chi_sq_min = min(chi_sq_min, chi_sq)
!  chi_sq_max = max(chi_sq_max, chi_sq)
!  SOFT_TOUCH chi_sq_min chi_sq_max
!END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, psi_norm_lr ]

  BEGIN_DOC
  ! < phi_l / phi_r > should be sampled with phi_r^2
  END_DOC

  implicit none

  psi_norm_lr = psidet_left_value / psidet_right_value

  psi_norm_lr_min = min(psi_norm_lr_min, psi_norm_lr)
  psi_norm_lr_max = max(psi_norm_lr_max, psi_norm_lr)
  SOFT_TOUCH psi_norm_lr_min psi_norm_lr_max

END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, sgn_lr ]

  BEGIN_DOC
  ! sign of Phi_L x Phi_R
  END_DOC

  implicit none
  double precision :: tmp

  tmp = psidet_left_value * psidet_right_value

  sgn_lr = tmp / dabs(tmp)

  sgn_lr_min = min(sgn_lr_min, sgn_lr)
  sgn_lr_max = max(sgn_lr_max, sgn_lr)
  SOFT_TOUCH sgn_lr_min sgn_lr_max

END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vartc_H ]

  BEGIN_DOC
  !
  ! < phi e^-J | H^2 | phi e^+J > = Eloc_mJpsi x Eloc_Jpsi
  !
  !  with :
  ! Eloc_Jpsi  = H (\Phi e^{+J}) / (\Phi e^{+J})
  ! Eloc_mJpsi = H (\Phi e^{-J}) / (\Phi e^{-J})
  !
  ! should be sampled with phi^2
  END_DOC

  implicit none

  vartc_H = Eloc_mJpsi * Eloc_Jpsi

  vartc_H_min = min(vartc_H_min, vartc_H)
  vartc_H_max = max(vartc_H_max, vartc_H)
  SOFT_TOUCH vartc_H_min vartc_H_max

END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

