
! ---

BEGIN_PROVIDER [ double precision, ci_Osvd, (size_ci_Osvd) ]

  BEGIN_DOC
  !
  ! Dimensions : n_svd_alpha * n_svd_beta * n_svd_alpha * n_svd_beta 
  END_DOC

  implicit none
  integer          :: k, kp, l, lp
  integer          :: ii0, ii1, ii2, ii3
  integer          :: iii, jjj
  double precision :: f

  iii = n_svd_alpha * n_svd_beta * n_svd_beta
  jjj = n_svd_alpha * n_svd_beta

  do k = 1, n_svd_alpha
    ii0 = (k-1) * iii 

    do kp = 1, n_svd_beta
      ii1 = ii0 + (kp-1) * jjj 

      f = det_alpha_value_SVD(k) * det_beta_value_SVD(kp) * psidet_right_inv

      do l = 1, n_svd_alpha
        ii2 = ii1 + (l-1) * n_svd_beta

        do lp = 1, n_svd_beta
          ii3 = ii2 + lp

          ci_Osvd(ii3) = det_alpha_value_SVD(l) * det_beta_value_SVD(lp) * psidet_right_inv * f
        enddo
      enddo
    enddo
  enddo

  ci_Osvd_min = min(ci_Osvd_min, minval(ci_Osvd))
  ci_Osvd_max = max(ci_Osvd_max, maxval(ci_Osvd))
  SOFT_TOUCH ci_Osvd_min ci_Osvd_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ci_Osvd_diag, (size_ci_Osvd_diag) ]

  BEGIN_DOC
  !
  ! Dimensions : n_svd_coefs * n_svd_coefs
  END_DOC

  implicit none
  integer          :: k, l, kk
  double precision :: f

  do k = 1, n_svd_coefs

    kk = (k-1) * n_svd_coefs 
    f  = det_alpha_value_SVD(k) * det_beta_value_SVD(k) * psidet_right_inv

    do l = 1, n_svd_coefs

      ci_Osvd_diag(kk+l) = det_alpha_value_SVD(l) * det_beta_value_SVD(l) * psidet_right_inv * f
    enddo
  enddo

  ci_Osvd_diag_min = min(ci_Osvd_diag_min, minval(ci_Osvd_diag))
  ci_Osvd_diag_max = max(ci_Osvd_diag_max, maxval(ci_Osvd_diag))
  SOFT_TOUCH ci_Osvd_diag_min ci_Osvd_diag_max
END_PROVIDER

! ---

