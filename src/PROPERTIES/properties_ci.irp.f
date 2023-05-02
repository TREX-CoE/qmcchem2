
! ---

BEGIN_PROVIDER [ double precision, ci_overlap_psidet, (size_ci_overlap_psidet) ]
  implicit none
  BEGIN_DOC
  ! < Phi_0 | det(j) >
  !
  ! Dimensions : det_num
  END_DOC

  integer                        :: i, j, k

  do k=1,det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    ci_overlap_psidet(k) = det_right_alpha_value(i)*det_right_beta_value (j)*psidet_right_inv
  enddo

  ci_overlap_psidet_min = min(ci_overlap_psidet_min,minval(ci_overlap_psidet))
  ci_overlap_psidet_max = max(ci_overlap_psidet_max,maxval(ci_overlap_psidet))
  SOFT_TOUCH ci_overlap_psidet_min ci_overlap_psidet_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ci_h_psidet, (size_ci_h_psidet) ]
  implicit none
  BEGIN_DOC
  ! < Phi_0 | H | det(j) >
  !
  ! Dimensions : det_num
  END_DOC

  integer                        :: i, j, k, l
  double precision :: T, tmp

  do k=1,det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    T = 0.d0
    do l=1,elec_alpha_num
      T += det_right_alpha_grad_lapl(4,l,i)*det_right_beta_value (j)
    enddo
    do l=elec_alpha_num+1,elec_num
      T += det_right_beta_grad_lapl (4,l,j)*det_right_alpha_value(i)
    enddo
    ci_h_psidet(k) = -0.5d0*T + E_pot * det_right_alpha_value(i)*det_right_beta_value (j)
    ci_h_psidet(k) *= psidet_right_inv
  enddo

  ci_h_psidet_min = min(ci_h_psidet_min,minval(ci_h_psidet))
  ci_h_psidet_max = max(ci_h_psidet_max,maxval(ci_h_psidet))
  SOFT_TOUCH ci_h_psidet_min ci_h_psidet_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ci_overlap_matrix, (size_ci_overlap_matrix) ]

  BEGIN_DOC
  ! < det(i) | det(j) >
  !
  ! Dimensions : det_num*det_num
  END_DOC

  implicit none
  integer          :: i, j, k, l, m, n, kk
  double precision :: f

  do k = 1, det_num

    kk = det_num * (k-1) 
    i  = det_coef_matrix_rows(k)
    j  = det_coef_matrix_columns(k)
    f  = det_right_alpha_value(i) * det_right_beta_value(j) * psidet_right_inv * psidet_right_inv

    do l = 1, det_num
      m = det_coef_matrix_rows(l)
      n = det_coef_matrix_columns(l)

      ci_overlap_matrix(kk+l) = det_right_alpha_value(m) * det_right_beta_value(n) * f
    enddo
  enddo

  ci_overlap_matrix_min = min(ci_overlap_matrix_min,minval(ci_overlap_matrix))
  ci_overlap_matrix_max = max(ci_overlap_matrix_max,maxval(ci_overlap_matrix))
  SOFT_TOUCH ci_overlap_matrix_min ci_overlap_matrix_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ci_h_matrix, (size_ci_h_matrix) ]

  BEGIN_DOC
  ! < det(i) |H| det(j) >
  !
  ! Dimensions : det_num*det_num
  END_DOC

  implicit none
  integer          :: i, j, k, l, ll, m, n, e, ee
  double precision :: f, g, h, T, V, j_lapl_inv

  ! (Lapl J)/J
  j_lapl_inv = 0.d0
  do e = 1, elec_num
    j_lapl_inv += jast_lapl_jast_inv(e)
  enddo

  do l = 1, det_num
    ll = det_num * (l-1)
    m  = det_coef_matrix_rows(l)
    n  = det_coef_matrix_columns(l)

    ! Lapl D 
    T = det_right_alpha_lapl_sum(m) * det_right_beta_value (n) &
      + det_right_beta_lapl_sum (n) * det_right_alpha_value(m) 

    if (j_lapl_inv /= 0.d0) then

      ! D (Lapl J)/J
      T += det_right_alpha_value(m) * det_right_beta_value(n) * j_lapl_inv

      ! 2 (grad D).(Grad J)/J
      g = 0.d0
      do e = 1, elec_alpha_num
        g += det_right_alpha_grad_lapl(1,e,m) * jast_grad_jast_inv_x(e) + &
             det_right_alpha_grad_lapl(2,e,m) * jast_grad_jast_inv_y(e) + &
             det_right_alpha_grad_lapl(3,e,m) * jast_grad_jast_inv_z(e)
      enddo
      h = 0.d0
      do e = 1, elec_beta_num
        ee = elec_alpha_num + e
        h += det_right_beta_grad_lapl(1,e,n) * jast_grad_jast_inv_x(ee) + &
             det_right_beta_grad_lapl(2,e,n) * jast_grad_jast_inv_y(ee) + &
             det_right_beta_grad_lapl(3,e,n) * jast_grad_jast_inv_z(ee)
      enddo
      T += 2.d0*( g * det_right_beta_value(n) + h * det_right_alpha_value(m) )

    endif

    g = det_right_alpha_value(m) * det_right_beta_value(n)
    V = (E_pot + E_nucl) * g
    if (do_pseudo) then
      do e = 1, elec_alpha_num
        V -= pseudo_right_non_local(e) * g
        V += det_right_alpha_pseudo(e,m) * det_right_beta_value(n)
      enddo
      do e = 1, elec_beta_num
        V -= pseudo_right_non_local(elec_alpha_num+e) * g
        !V -= pseudo_right_non_local(e) * g
        V += det_right_alpha_value(m) * det_right_beta_pseudo(e,n)
      enddo
    endif 

    f = -0.5d0 * T + V
    f *= psidet_right_inv * psidet_right_inv

    do k = 1, det_num
      i = det_coef_matrix_rows(k)
      j = det_coef_matrix_columns(k)
      ci_h_matrix(ll + k) = f * det_right_alpha_value(i) * det_right_beta_value(j)
    enddo
  enddo

  ci_h_matrix_min = min(ci_h_matrix_min, minval(ci_h_matrix))
  ci_h_matrix_max = max(ci_h_matrix_max, maxval(ci_h_matrix))
  SOFT_TOUCH ci_h_matrix_min ci_h_matrix_max
END_PROVIDER

! ---

!BEGIN_PROVIDER [ double precision, H_DJ, (size_H_DJ) ]
!
!  BEGIN_DOC
!  ! H | det e^J > / [ phi e^{J} ]
!  !
!  ! Dimensions : det_num
!  END_DOC
!
!  implicit none
!  integer          :: l, m, n, e, ee
!  double precision :: g, h, T, V, j_lapl_inv
!
!  ! (Lapl J)/J
!  j_lapl_inv = 0.d0
!  do e = 1, elec_num
!    j_lapl_inv += jast_lapl_jast_inv(e)
!  enddo
!
!  do l = 1, det_num
!    m  = det_coef_matrix_rows(l)
!    n  = det_coef_matrix_columns(l)
!
!    ! Lapl D 
!    T = det_right_alpha_lapl_sum(m) * det_right_beta_value (n) &
!      + det_right_beta_lapl_sum (n) * det_right_alpha_value(m) 
!
!    if (j_lapl_inv /= 0.d0) then
!
!      ! D (Lapl J)/J
!      T += det_right_alpha_value(m) * det_right_beta_value(n) * j_lapl_inv
!
!      ! 2 (grad D).(Grad J)/J
!      g = 0.d0
!      do e = 1, elec_alpha_num
!        g += det_right_alpha_grad_lapl(1,e,m) * jast_grad_jast_inv_x(e) + &
!             det_right_alpha_grad_lapl(2,e,m) * jast_grad_jast_inv_y(e) + &
!             det_right_alpha_grad_lapl(3,e,m) * jast_grad_jast_inv_z(e)
!      enddo
!      h = 0.d0
!      do e = 1, elec_beta_num
!        ee = elec_alpha_num + e
!        h += det_right_beta_grad_lapl(1,e,n) * jast_grad_jast_inv_x(ee) + &
!             det_right_beta_grad_lapl(2,e,n) * jast_grad_jast_inv_y(ee) + &
!             det_right_beta_grad_lapl(3,e,n) * jast_grad_jast_inv_z(ee)
!      enddo
!      T += 2.d0*( g * det_right_beta_value(n) + h * det_right_alpha_value(m) )
!
!    endif
!
!    g = det_right_alpha_value(m) * det_right_beta_value(n)
!    V = (E_pot + E_nucl) * g
!    if (do_pseudo) then
!      do e = 1, elec_alpha_num
!        V -= pseudo_right_non_local(e) * g
!        V += det_right_alpha_pseudo(e,m) * det_right_beta_value(n)
!      enddo
!      do e = 1, elec_beta_num
!        V -= pseudo_right_non_local(elec_alpha_num+e) * g
!        V += det_right_alpha_value(m) * det_right_beta_pseudo(e,n)
!      enddo
!    endif 
!
!    H_DJ(l) = ( -0.5d0*T + V ) * psidet_right_inv
!  enddo
!
!  H_DJ_min = min(H_DJ_min, minval(H_DJ))
!  H_DJ_max = max(H_DJ_max, maxval(H_DJ))
!  SOFT_TOUCH H_DJ_min H_DJ_max
!END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ci_h_matrix_diag, (size_ci_h_matrix_diag) ]
  implicit none
  BEGIN_DOC
  ! < det(i) |H| det(j) >
  !
  ! Dimensions : det_num
  END_DOC

  integer                        :: i, j, k, l, m, n, e
  double precision :: f, g, h, T, V

  do l=1,det_num
      m = det_coef_matrix_rows(l)
      n = det_coef_matrix_columns(l)
      ! Lapl D 
      g = 0.d0
      do e=1,elec_alpha_num
        g += det_right_alpha_grad_lapl(4,e,m) * det_right_beta_value (n)
      enddo
      do e=elec_alpha_num+1,elec_num
        g += det_right_alpha_value(m) * det_right_beta_grad_lapl(4,e,n)
      enddo
      T = g
      ! D (Lapl J)/J
      g = 0.d0
      do e=1,elec_num
        g += jast_lapl_jast_inv(e)
      enddo
      T += det_right_alpha_value(m) * det_right_beta_value(n) * g
      ! 2 (grad D).(Grad J)/J
      g = 0.d0
      do e=1,elec_alpha_num
        g += &
          det_right_alpha_grad_lapl(1,e,m) * jast_grad_jast_inv_x(e) + &
          det_right_alpha_grad_lapl(2,e,m) * jast_grad_jast_inv_y(e) + &
          det_right_alpha_grad_lapl(3,e,m) * jast_grad_jast_inv_z(e)
      enddo
      h = 0.d0
      do e=elec_alpha_num+1,elec_num
        h += &
          det_right_beta_grad_lapl(1,e,n) * jast_grad_jast_inv_x(e) + &
          det_right_beta_grad_lapl(2,e,n) * jast_grad_jast_inv_y(e) + &
          det_right_beta_grad_lapl(3,e,n) * jast_grad_jast_inv_z(e)
      enddo
      T += 2.d0*( g * det_right_beta_value(n) + h * det_right_alpha_value(m) )
      g = det_right_alpha_value(m)*det_right_beta_value(n)
      V = E_pot* g
      do e=1,elec_alpha_num
        V -= pseudo_right_non_local(e)* g
        V += det_right_alpha_pseudo(e,m) * det_right_beta_value(n)
      enddo
      do e=elec_alpha_num+1,elec_num
        V -= pseudo_right_non_local(e)* g
        V += det_right_alpha_value(m) * det_right_beta_pseudo(e,n)
      enddo
      f = -0.5d0*T + V
      f *= psidet_right_inv * psidet_right_inv
      ci_h_matrix_diag(l) = f * det_right_alpha_value(m)*det_right_beta_value (n)
  enddo

  ci_h_matrix_diag_min = min(ci_h_matrix_diag_min,minval(ci_h_matrix_diag))
  ci_h_matrix_diag_max = max(ci_h_matrix_diag_max,maxval(ci_h_matrix_diag))
  SOFT_TOUCH ci_h_matrix_diag_min ci_h_matrix_diag_max
END_PROVIDER

! ---


