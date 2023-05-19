program SVD_info 

  PROVIDE ezfio_filename

  implicit none

  print *, ' Number of determinants                   : ', det_num
  print *, ' Number of unique alpha/beta determinants : ', det_alpha_num, det_beta_num
 
  print*, ' SVD param:'
  print*, ' n_svd_alpha = ', n_svd_alpha
  print*, ' n_svd_beta = ', n_svd_beta
  print*, ' la  = ', la 
  print*, ' lb  = ', lb 
  print*, ' ld  = ', ld 
  print*, ' lla = ', lla 
  print*, ' llb = ', llb 
  print*, ' lld = ', lld 
  print*, ' n_dress_svd = ', n_dress_svd

  call t_info() 

  !call check_UV_norm()

end

! ---

subroutine check_UV_norm()

  implicit none
  integer                       :: i, j
  double precision              :: norm_mat 
  double precision, allocatable :: UUt(:,:), VVt(:,:)

  ! ---

  print*, ' computing U x U.t ... '

  allocate( UUt(det_alpha_num,det_alpha_num) )
  call dgemm( 'N', 'T', det_alpha_num, det_alpha_num, n_svd_alpha, 1.d0 &
             , psi_svd_alpha(1,1,1), size(psi_svd_alpha, 1)             &
             , psi_svd_alpha(1,1,1), size(psi_svd_alpha, 1)             &
             , 0.d0, UUt, size(UUt, 1) )

  print*, ' norm U x U.t : '

  norm_mat = 0.d0
  do i = 1, det_alpha_num
    norm_mat = norm_mat + UUt(i,i) * UUt(i,i)
  enddo
  print*, ' diag elements = ', dsqrt(norm_mat)

  norm_mat = 0.d0
  do i = 1, det_alpha_num
    do j = 1, det_alpha_num
      if( i.eq.j ) cycle
      norm_mat = norm_mat + UUt(j,i) * UUt(j,i)
    enddo
  enddo
  print*, ' extra-diag elements = ', dsqrt(norm_mat)

  deallocate( UUt )

  ! ---

  print*, ' computing V x V.t ... '

  allocate( VVt(det_beta_num,det_beta_num) )
  call dgemm( 'N', 'T', det_beta_num, det_beta_num, n_svd_beta, 1.d0 &
             , psi_svd_beta(1,1,1), size(psi_svd_beta, 1)            &
             , psi_svd_beta(1,1,1), size(psi_svd_beta, 1)            &
             , 0.d0, VVt, size(VVt, 1) )

  print*, ' norm V x V.t : '
  norm_mat = 0.d0
  do i = 1, det_beta_num
    norm_mat = norm_mat + VVt(i,i) * VVt(i,i)
  enddo
  print*, ' diag elements = ', dsqrt(norm_mat)
  norm_mat = 0.d0
  do i = 1, det_beta_num
    do j = 1, det_beta_num
      if( i.eq.j ) cycle
      norm_mat = norm_mat + VVt(j,i) * VVt(j,i)
    enddo
  enddo
  print*, ' extra-diag elements = ', dsqrt(norm_mat)

  deallocate( VVt )

  ! ---

end subroutine check_UV_norm

! ---

subroutine t_info()

  implicit none
  integer          :: i, imax
  double precision :: ti, tf 

  call CPU_TIME(ti)
  PROVIDE psi_svd_alpha psi_svd_beta
  call CPU_TIME(tf)
  print*, 'read SVD matrices after (sec)', tf-ti

  PROVIDE E_loc
  print *,  'E_loc = ', E_loc
  print *, 'computing det_alpha_value_svd & det_beta_value_svd ...'
  PROVIDE det_alpha_value_svd det_beta_value_svd 

  call CPU_TIME(ti)
  imax = 100
  do i = 1, imax
    PROVIDE E_loc
    PROVIDE det_alpha_value_svd det_beta_value_svd 
    TOUCH elec_coord
  enddo
  call CPU_TIME(tf)
  print *, 'Time for the calculation of E_loc, det_alpha_value_svd & det_beta_value_svd (sec) = ', (tf-ti)/float(imax)

end subroutine t_info

! ---
