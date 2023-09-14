! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJsimple_m_Hmu_no3b_svd, (size_vec_HJsimple_m_Hmu_no3b_svd) ]

  BEGIN_DOC
  !
  ! Dimensions : n_dress_svd 
  END_DOC

  implicit none
  integer                       :: ic, i, j, k
  double precision              :: w
  double precision, allocatable :: tmp_beta(:)

  allocate( tmp_beta(det_beta_num) )

  w = EJsimple_m_EJmu1b_no3b * Psi_value_inv * jast_value_inv
  tmp_beta(:) = w * det_beta_value_svd(:) 

  ic = 0
  do i = 1, ld 
    w = det_alpha_value_svd(i)
    do j = 1, ld
      vec_HJsimple_m_Hmu_no3b_svd(ic+j) = w * tmp_beta(j) 
    enddo
    ic = ic+ld
  enddo

  do i = ld+1, lld 
    vec_HJsimple_m_Hmu_no3b_svd(ic+i-ld) = det_alpha_value_svd(i) * tmp_beta(i) 
  enddo
  ic = ic + lld - ld

  do i = ld+1, lla 
    w = det_alpha_value_svd(i) 
    do j = 1, lb
      vec_HJsimple_m_Hmu_no3b_svd(ic+j) = w * tmp_beta(j) 
    enddo
    ic = ic+lb
  enddo

  do i = 1, la 
    w = det_alpha_value_svd(i) 
    do j = ld+1, llb
      vec_HJsimple_m_Hmu_no3b_svd(ic+j-ld) = w * tmp_beta(j) 
    enddo
    ic = ic+llb-ld
  enddo

  deallocate( tmp_beta )

  vec_HJsimple_m_Hmu_no3b_svd_min = min(vec_HJsimple_m_Hmu_no3b_svd_min, minval(vec_HJsimple_m_Hmu_no3b_svd))
  vec_HJsimple_m_Hmu_no3b_svd_max = max(vec_HJsimple_m_Hmu_no3b_svd_max, maxval(vec_HJsimple_m_Hmu_no3b_svd))
  SOFT_TOUCH vec_HJsimple_m_Hmu_no3b_svd_min vec_HJsimple_m_Hmu_no3b_svd_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJsimple_m_H_svd, (size_vec_HJsimple_m_H_svd) ]

  BEGIN_DOC
  !
  ! Dimensions : n_dress_svd 
  END_DOC

  implicit none
  integer                       :: ic, i, j, k
  double precision              :: w
  double precision, allocatable :: tmp_beta(:)

  allocate( tmp_beta(det_beta_num) )

  w = deltaE_Jsimple * Psi_value_inv * jast_value_inv
  tmp_beta(:) = w * det_beta_value_svd(:) 

  ic = 0
  do i = 1, ld 
    w = det_alpha_value_svd(i)
    do j = 1, ld
      vec_HJsimple_m_H_svd(ic+j) = w * tmp_beta(j) 
    enddo
    ic = ic+ld
  enddo

  do i = ld+1, lld 
    vec_HJsimple_m_H_svd(ic+i-ld) = det_alpha_value_svd(i) * tmp_beta(i) 
  enddo
  ic = ic + lld - ld

  do i = ld+1, lla 
    w = det_alpha_value_svd(i) 
    do j = 1, lb
      vec_HJsimple_m_H_svd(ic+j) = w * tmp_beta(j) 
    enddo
    ic = ic+lb
  enddo

  do i = 1, la 
    w = det_alpha_value_svd(i) 
    do j = ld+1, llb
      vec_HJsimple_m_H_svd(ic+j-ld) = w * tmp_beta(j) 
    enddo
    ic = ic+llb-ld
  enddo

  deallocate( tmp_beta )

  vec_HJsimple_m_H_svd_min = min(vec_HJsimple_m_H_svd_min, minval(vec_HJsimple_m_H_svd))
  vec_HJsimple_m_H_svd_max = max(vec_HJsimple_m_H_svd_max, maxval(vec_HJsimple_m_H_svd))
  SOFT_TOUCH vec_HJsimple_m_H_svd_min vec_HJsimple_m_H_svd_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

! -------------------------------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, vec_HJmu_m_H_svd, (size_vec_HJmu_m_H_svd) ]

  BEGIN_DOC
  !
  ! Dimensions : n_dress_svd 
  END_DOC

  implicit none
  integer                       :: ic, i, j, k
  double precision              :: w
  double precision, allocatable :: tmp_beta(:)

  allocate( tmp_beta(det_beta_num) )

  w = deltaE_Jmu1b * Psi_value_inv * jast_value_inv
  tmp_beta(:) = w * det_beta_value_svd(:) 

  ic = 0
  do i = 1, ld 
    w = det_alpha_value_svd(i)
    do j = 1, ld
      vec_HJmu_m_H_svd(ic+j) = w * tmp_beta(j) 
    enddo
    ic = ic+ld
  enddo

  do i = ld+1, lld 
    vec_HJmu_m_H_svd(ic+i-ld) = det_alpha_value_svd(i) * tmp_beta(i) 
  enddo
  ic = ic + lld - ld

  do i = ld+1, lla 
    w = det_alpha_value_svd(i) 
    do j = 1, lb
      vec_HJmu_m_H_svd(ic+j) = w * tmp_beta(j) 
    enddo
    ic = ic+lb
  enddo

  do i = 1, la 
    w = det_alpha_value_svd(i) 
    do j = ld+1, llb
      vec_HJmu_m_H_svd(ic+j-ld) = w * tmp_beta(j) 
    enddo
    ic = ic+llb-ld
  enddo

  deallocate( tmp_beta )

  vec_HJmu_m_H_svd_min = min(vec_HJmu_m_H_svd_min, minval(vec_HJmu_m_H_svd))
  vec_HJmu_m_H_svd_max = max(vec_HJmu_m_H_svd_max, maxval(vec_HJmu_m_H_svd))
  SOFT_TOUCH vec_HJmu_m_H_svd_min vec_HJmu_m_H_svd_max
END_PROVIDER

! -------------------------------------------------------------------------------------------------

! ---

