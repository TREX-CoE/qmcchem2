
! ---

BEGIN_PROVIDER [ double precision, DJ_Eloc_DJ, (size_DJ_Eloc_DJ) ]

  BEGIN_DOC
  !
  ! DJ_Eloc_DJ = < D(i) e^{J} | E_loc | D(j) e^{J} >
  !            = < [D(i)/Phi] * [D(j)/Phi] * E_loc >_{Psi^2}
  !
  ! Dimensions : det_num*det_num
  END_DOC

  implicit none
  integer          :: i, j, k, l, m, n, ll
  double precision :: wl, wk

  do l = 1, det_num
    ll = det_num * (l-1)
    m  = det_coef_matrix_rows(l)
    n  = det_coef_matrix_columns(l)
    wl = det_right_alpha_value(m) * det_right_beta_value(n) * psidet_right_inv
    do k = 1, det_num
      i  = det_coef_matrix_rows(k)
      j  = det_coef_matrix_columns(k)
      wk = det_right_alpha_value(i) * det_right_beta_value(j) * psidet_right_inv

      DJ_Eloc_DJ(ll) = E_loc * wl * wk
    enddo
  enddo

  DJ_Eloc_DJ_min = min(DJ_Eloc_DJ_min, minval(DJ_Eloc_DJ))
  DJ_Eloc_DJ_max = max(DJ_Eloc_DJ_max, maxval(DJ_Eloc_DJ))
  SOFT_TOUCH DJ_Eloc_DJ_min DJ_Eloc_DJ_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, DJ_Eloc_DJ_diag, (size_DJ_Eloc_DJ_diag) ]

  BEGIN_DOC
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: l, m, n
  double precision :: wl

  do l = 1, det_num
    m  = det_coef_matrix_rows(l)
    n  = det_coef_matrix_columns(l)
    wl = det_right_alpha_value(m) * det_right_beta_value(n) * psidet_right_inv

    DJ_Eloc_DJ_diag(l) = E_loc * wl * wl
  enddo

  DJ_Eloc_DJ_diag_min = min(DJ_Eloc_DJ_diag_min, minval(DJ_Eloc_DJ_diag))
  DJ_Eloc_DJ_diag_max = max(DJ_Eloc_DJ_diag_max, maxval(DJ_Eloc_DJ_diag))
  SOFT_TOUCH DJ_Eloc_DJ_diag_min DJ_Eloc_DJ_diag_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, DJ_Jphi, (size_DJ_Jphi) ]

  BEGIN_DOC
  !
  ! DJ_Jphi(i) = < D(i) e^{J} | e^{J} Phi >
  !            = < D(i) / Phi >_{Psi^2}
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k

  do k = 1, det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    DJ_Jphi(k) = psidet_right_inv * det_right_alpha_value(i) * det_right_beta_value(j)
  enddo

  DJ_Jphi_min = min(DJ_Jphi_min, minval(DJ_Jphi))
  DJ_Jphi_max = max(DJ_Jphi_max, maxval(DJ_Jphi))
  SOFT_TOUCH DJ_Jphi_min DJ_Jphi_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, DJ_H_Jphi, (size_DJ_H_Jphi) ]

  BEGIN_DOC
  !
  ! DJ_H_Jphi(i) = < D(i) e^{J} | H | e^{J} Phi >
  !              = < D(i) E_loc / Phi >_{Psi^2}
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = E_loc * psidet_right_inv
  do k = 1, det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    DJ_H_Jphi(k) = w * det_right_alpha_value(i) * det_right_beta_value(j)
  enddo

  DJ_H_Jphi_min = min(DJ_H_Jphi_min, minval(DJ_H_Jphi))
  DJ_H_Jphi_max = max(DJ_H_Jphi_max, maxval(DJ_H_Jphi))
  SOFT_TOUCH DJ_H_Jphi_min DJ_H_Jphi_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, dcI_Eloc, (size_dcI_Eloc) ]

  BEGIN_DOC
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, k
  double precision :: w

  w = E_loc * psidet_right_inv
  do k = 1, det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)

    dcI_Eloc(k) = H_DJ(k) - w * det_right_alpha_value(i) * det_right_beta_value(j)
  enddo

  dcI_Eloc_min = min(dcI_Eloc_min, minval(dcI_Eloc))
  dcI_Eloc_max = max(dcI_Eloc_max, maxval(dcI_Eloc))
  SOFT_TOUCH dcI_Eloc_min dcI_Eloc_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, DJ_dcI_Eloc_DJ, (size_DJ_dcI_Eloc_DJ) ]

  BEGIN_DOC
  !
  ! DJ_dcI_Eloc_DJ(i,j) = < D(i) e^{J} | \partial_{C_I} Eloc >
  !                     = < D(i) dcI_Eloc(j) / Phi >_{Psi^2}
  !
  ! Dimensions : det_num*det_num
  END_DOC

  implicit none
  integer          :: i, j, k, l, ll
  double precision :: w

  do l = 1, det_num
    ll = det_num * (l-1)
    i = det_coef_matrix_rows(l)
    j = det_coef_matrix_columns(l)
    w = psidet_right_inv * det_right_alpha_value(i) * det_right_beta_value(j)
    do k = 1, det_num
      DJ_dcI_Eloc_DJ(ll + k) = w * dcI_Eloc(k)
    enddo
  enddo

  DJ_dcI_Eloc_DJ_min = min(DJ_dcI_Eloc_DJ_min, minval(DJ_dcI_Eloc_DJ))
  DJ_dcI_Eloc_DJ_max = max(DJ_dcI_Eloc_DJ_max, maxval(DJ_dcI_Eloc_DJ))
  SOFT_TOUCH DJ_dcI_Eloc_DJ_min DJ_dcI_Eloc_DJ_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, DJ_dcI_Eloc_DJ_diag, (size_DJ_dcI_Eloc_DJ_diag) ]

  BEGIN_DOC
  !
  ! Dimensions : det_num
  END_DOC

  implicit none
  integer          :: i, j, l
  double precision :: w

  do l = 1, det_num
    i = det_coef_matrix_rows(l)
    j = det_coef_matrix_columns(l)
    w = psidet_right_inv * det_right_alpha_value(i) * det_right_beta_value(j)

    DJ_dcI_Eloc_DJ_diag(l) = w * dcI_Eloc(l)
  enddo

  DJ_dcI_Eloc_DJ_diag_min = min(DJ_dcI_Eloc_DJ_diag_min, minval(DJ_dcI_Eloc_DJ_diag))
  DJ_dcI_Eloc_DJ_diag_max = max(DJ_dcI_Eloc_DJ_diag_max, maxval(DJ_dcI_Eloc_DJ_diag))
  SOFT_TOUCH DJ_dcI_Eloc_DJ_diag_min DJ_dcI_Eloc_DJ_diag_max
END_PROVIDER

! ---


