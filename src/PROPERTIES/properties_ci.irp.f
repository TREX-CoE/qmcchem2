BEGIN_PROVIDER [ double precision, psi_norm ]
 implicit none
  BEGIN_DOC
  ! <1/J^2>
  END_DOC
  psi_norm = jast_value_inv*jast_value_inv

  psi_norm_min = min(psi_norm_min,psi_norm)
  psi_norm_max = max(psi_norm_max,psi_norm)
  SOFT_TOUCH psi_norm_min psi_norm_max
END_PROVIDER

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
    ci_overlap_psidet(k) = det_alpha_value(i)*det_beta_value (j)*psidet_inv
  enddo

  ci_overlap_psidet_min = min(ci_overlap_psidet_min,minval(ci_overlap_psidet))
  ci_overlap_psidet_max = max(ci_overlap_psidet_max,maxval(ci_overlap_psidet))
  SOFT_TOUCH ci_overlap_psidet_min ci_overlap_psidet_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, ci_h_psidet, (size_ci_h_psidet) ]
  implicit none
  BEGIN_DOC
  ! < Phi_0 | H | det(j) >
  !
  ! Dimensions : det_num
  END_DOC

  integer                        :: i, j, k, l
  double precision :: T

  do k=1,det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    T = det_alpha_lapl_sum(i)*det_beta_value(j) + det_beta_lapl_sum(j)*det_alpha_value(i)
    ci_h_psidet(k) = -0.5d0*T + (E_pot + E_nucl) * det_alpha_value(i)*det_beta_value (j)
    ci_h_psidet(k) *= psi_value_inv * jast_value_inv
  enddo

  ci_h_psidet_min = min(ci_h_psidet_min,minval(ci_h_psidet))
  ci_h_psidet_max = max(ci_h_psidet_max,maxval(ci_h_psidet))
  SOFT_TOUCH ci_h_psidet_min ci_h_psidet_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, ci_overlap_matrix, (size_ci_overlap_matrix) ]
  implicit none
  BEGIN_DOC
  ! < det(i) | det(j) >
  !
  ! Dimensions : det_num*det_num
  END_DOC

  integer                        :: i, j, k, l, m, n
  double precision :: f

  do k=1,det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    f = det_alpha_value(i)*det_beta_value (j)*psidet_inv*psidet_inv
    do l=1,det_num
      m = det_coef_matrix_rows(l)
      n = det_coef_matrix_columns(l)
      ci_overlap_matrix( det_num*(k-1) + l) = det_alpha_value(m)*det_beta_value(n) * f
    enddo
  enddo

  ci_overlap_matrix_min = min(ci_overlap_matrix_min,minval(ci_overlap_matrix))
  ci_overlap_matrix_max = max(ci_overlap_matrix_max,maxval(ci_overlap_matrix))
  SOFT_TOUCH ci_overlap_matrix_min ci_overlap_matrix_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, ci_h_matrix, (size_ci_h_matrix) ]
  implicit none
  BEGIN_DOC
  ! < det(i) |H| det(j) >
  !
  ! Dimensions : det_num*det_num
  END_DOC

  integer                        :: i, j, k, l, m, n, e
  double precision :: f, g, h, T, V, j_lapl_inv

  ! (Lapl J)/J
  j_lapl_inv = 0.d0
  do e=1,elec_num
    j_lapl_inv += jast_lapl_jast_inv(e)
  enddo

  do l=1,det_num
      m = det_coef_matrix_rows(l)
      n = det_coef_matrix_columns(l)
      ! Lapl D 
      T = det_alpha_lapl_sum(m) * det_beta_value (n) &
          + det_alpha_value(m) * det_beta_lapl_sum(n)
      if (j_lapl_inv /= 0.d0) then
        ! D (Lapl J)/J
        T += det_alpha_value(m) * det_beta_value(n) * j_lapl_inv

        ! 2 (grad D).(Grad J)/J
        g = 0.d0
        do e=1,elec_alpha_num
          g += &
            det_alpha_grad_lapl(1,e,m) * jast_grad_jast_inv_x(e) + &
            det_alpha_grad_lapl(2,e,m) * jast_grad_jast_inv_y(e) + &
            det_alpha_grad_lapl(3,e,m) * jast_grad_jast_inv_z(e)
        enddo
        h = 0.d0
        do e=1,elec_beta_num
          h += &
            det_beta_grad_lapl(1,e,n) * jast_grad_jast_inv_x(elec_alpha_num+e) + &
            det_beta_grad_lapl(2,e,n) * jast_grad_jast_inv_y(elec_alpha_num+e) + &
            det_beta_grad_lapl(3,e,n) * jast_grad_jast_inv_z(elec_alpha_num+e)
        enddo
        T += 2.d0*( g * det_beta_value(n) + h * det_alpha_value(m) )
      endif
      g = det_alpha_value(m)*det_beta_value(n)
      V = (E_pot + E_nucl)* g
      if (do_pseudo) then
        do e=1,elec_alpha_num
          V -= pseudo_non_local(e)* g
          V += det_alpha_pseudo(e,m) * det_beta_value(n)
        enddo
        do e=1,elec_beta_num
          V -= pseudo_non_local(e)* g
          V += det_alpha_value(m) * det_beta_pseudo(e,n)
        enddo
      endif 
      f = -0.5d0*T + V
      f *= psidet_inv * psidet_inv
      do k=1,det_num
        i = det_coef_matrix_rows(k)
        j = det_coef_matrix_columns(k)
        ci_h_matrix( det_num*(l-1) + k) = f * &
           det_alpha_value(i)*det_beta_value (j)
      enddo
  enddo
  ci_h_matrix_min = min(ci_h_matrix_min,minval(ci_h_matrix))
  ci_h_matrix_max = max(ci_h_matrix_max,maxval(ci_h_matrix))
  SOFT_TOUCH ci_h_matrix_min ci_h_matrix_max
END_PROVIDER


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
        g += det_alpha_grad_lapl(4,e,m) * det_beta_value (n)
      enddo
      do e=1,elec_beta_num
        g += det_alpha_value(m) * det_beta_grad_lapl(4,e,n)
      enddo
      T = g
      ! D (Lapl J)/J
      g = 0.d0
      do e=1,elec_num
        g += jast_lapl_jast_inv(e)
      enddo
      T += det_alpha_value(m) * det_beta_value(n) * g
      ! 2 (grad D).(Grad J)/J
      g = 0.d0
      do e=1,elec_alpha_num
        g += &
          det_alpha_grad_lapl(1,e,m) * jast_grad_jast_inv_x(e) + &
          det_alpha_grad_lapl(2,e,m) * jast_grad_jast_inv_y(e) + &
          det_alpha_grad_lapl(3,e,m) * jast_grad_jast_inv_z(e)
      enddo
      h = 0.d0
      do e=1,elec_beta_num
        h += &
          det_beta_grad_lapl(1,e,n) * jast_grad_jast_inv_x(elec_alpha_num+e) + &
          det_beta_grad_lapl(2,e,n) * jast_grad_jast_inv_y(elec_alpha_num+e) + &
          det_beta_grad_lapl(3,e,n) * jast_grad_jast_inv_z(elec_alpha_num+e)
      enddo
      T += 2.d0*( g * det_beta_value(n) + h * det_alpha_value(m) )
      g = det_alpha_value(m)*det_beta_value(n)
      V = E_pot* g
      if (do_pseudo) then
        do e=1,elec_alpha_num
          V -= pseudo_non_local(e)* g
          V += det_alpha_pseudo(e,m) * det_beta_value(n)
        enddo
        do e=1,elec_beta_num
          V -= pseudo_non_local(e)* g
          V += det_alpha_value(m) * det_beta_pseudo(e,n)
        enddo
      endif
      f = -0.5d0*T + V
      f *= psidet_inv * psidet_inv
      ci_h_matrix_diag(l) = f * &
           det_alpha_value(m)*det_beta_value (n)
  enddo

  ci_h_matrix_diag_min = min(ci_h_matrix_diag_min,minval(ci_h_matrix_diag))
  ci_h_matrix_diag_max = max(ci_h_matrix_diag_max,maxval(ci_h_matrix_diag))
  SOFT_TOUCH ci_h_matrix_diag_min ci_h_matrix_diag_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, ci_h_transcor_psi, (size_ci_h_transcor_psi) ]
  implicit none
  BEGIN_DOC
  ! < det(i) e^{-J} |H| Psi >
  !
  ! Dimensions : det_num
  END_DOC

  integer                        :: i, j, k

  do k=1,det_num
    i = det_coef_matrix_rows(k)
    j = det_coef_matrix_columns(k)
    ci_h_transcor_psi(k) = E_loc * jast_value_inv * &
       det_alpha_value(i)*det_beta_value(j) * psi_value_inv
  enddo

  ci_h_transcor_psi_min = min(ci_h_transcor_psi_min,minval(ci_h_transcor_psi))
  ci_h_transcor_psi_max = max(ci_h_transcor_psi_max,maxval(ci_h_transcor_psi))
  SOFT_TOUCH ci_h_transcor_psi_min ci_h_transcor_psi_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, ci_dress, (size_ci_dress) ]
  implicit none
  BEGIN_DOC
  ! < det(i) e^{-J} |H| Psi >
  !
  ! Dimensions : det_num
  END_DOC

  integer          :: i, j, k, l
  double precision :: T, h_psidet, dij, f, E_noJ, dE

  h_psidet = -0.5d0*psidet_lapl*psidet_inv + E_pot + E_nucl 
  E_noJ = h_psidet 
  dE = E_loc - E_noJ
  do k=1,det_num
     i = det_coef_matrix_rows(k)
     j = det_coef_matrix_columns(k)
     f = det_alpha_value(i)*det_beta_value(j) * psi_value_inv * jast_value_inv
     ci_dress(k) = dE * f
  enddo

return
  integer                        :: m, n, e
  double precision :: g, h, V, j_lapl_inv, det_ab

  ! (Lapl J)/J
  j_lapl_inv = 0.d0
  do e=1,elec_num
    j_lapl_inv += jast_lapl_jast_inv(e)
  enddo

  do l=1,det_num
      m = det_coef_matrix_rows(l)
      n = det_coef_matrix_columns(l)
      ! Lapl D 
!      T = det_alpha_lapl_sum(m) * det_beta_value (n) &
!          + det_alpha_value(m) * det_beta_lapl_sum(n)
!      det_ab = det_alpha_value(m)*det_beta_value(n)
!      ci_dress(l) = -0.5d0*T + (E_pot + E_nucl) * det_ab
      T = 0.d0
      ci_dress(l) = 0.d0

      ! D (Lapl J)/J
      T += det_alpha_value(m) * det_beta_value(n) * j_lapl_inv

      ! 2 (grad D).(Grad J)/J
      g = 0.d0
      do e=1,elec_alpha_num
        g += &
          det_alpha_grad_lapl(1,e,m) * jast_grad_jast_inv_x(e) + &
          det_alpha_grad_lapl(2,e,m) * jast_grad_jast_inv_y(e) + &
          det_alpha_grad_lapl(3,e,m) * jast_grad_jast_inv_z(e)
      enddo
      h = 0.d0
      do e=1,elec_beta_num
        h += &
          det_beta_grad_lapl(1,e,n) * jast_grad_jast_inv_x(elec_alpha_num+e) + &
          det_beta_grad_lapl(2,e,n) * jast_grad_jast_inv_y(elec_alpha_num+e) + &
          det_beta_grad_lapl(3,e,n) * jast_grad_jast_inv_z(elec_alpha_num+e)
      enddo
      T += 2.d0*( g * det_beta_value(n) + h * det_alpha_value(m) )

      V = 0.d0 ! (E_pot + E_nucl)* det_ab
      if (do_pseudo) then
        do e=1,elec_alpha_num
          V -= pseudo_non_local(e)* det_ab
          V += det_alpha_pseudo(e,m) * det_beta_value(n)
        enddo
        do e=1,elec_beta_num
          V -= pseudo_non_local(e)* det_ab
          V += det_alpha_value(m) * det_beta_pseudo(e,n)
        enddo
      endif 
      f = -0.5d0*T + V !- ci_dress(l)
      ci_dress(l) = f *  psi_value_inv * jast_value_inv * jast_value_inv

  enddo


  ci_dress_min = min(ci_dress_min,minval(ci_dress))
  ci_dress_max = max(ci_dress_max,maxval(ci_dress))
  SOFT_TOUCH ci_dress_min ci_dress_max
END_PROVIDER

BEGIN_PROVIDER [ double precision, ci_dress_opt ]
  BEGIN_DOC
  ! Use for optimizing mu
  END_DOC
  implicit none
  integer          :: i, j, k, l
  double precision :: T, dij, f, E_noJ, dE
  ! energy = H \Phi / \Phi
  E_noJ = -0.5d0*psidet_lapl*psidet_inv + E_pot + E_nucl 
  dE = (E_loc - E_noJ) * psi_value_inv * jast_value_inv ! PsiJ.J
  k = 1
  i = det_coef_matrix_rows(   k)
  j = det_coef_matrix_columns(k)
  f = det_alpha_value(i) * det_beta_value(j)
  ci_dress_opt = dE * f
  ci_dress_opt_min = min(ci_dress_opt_min, ci_dress_opt)
  ci_dress_opt_max = max(ci_dress_opt_max, ci_dress_opt)
  SOFT_TOUCH ci_dress_opt_min ci_dress_opt_max
END_PROVIDER

BEGIN_PROVIDER [ double precision, ci_dress_Htilde, (size_ci_dress) ]
  implicit none
  BEGIN_DOC
  ! < det(i) e^{-J} |H| Psi >
  !
  ! Dimensions : det_num
  END_DOC

  integer          :: i, j, k, l
  double precision :: T, h_psidet, dij, f, E_noJ, dE

  E_noJ = -0.5d0*psidet_lapl*psidet_inv + E_pot + E_nucl 
  dE = E_loc - E_noJ
  do k=1,det_num
     i = det_coef_matrix_rows(k)
     j = det_coef_matrix_columns(k)
     f = det_alpha_value(i)*det_beta_value(j) * psi_value_inv * jast_value_inv
     ci_dress(k) = dE * f
  enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, ci_dress_H, (size_ci_dress) ]
  implicit none
  BEGIN_DOC
  ! < det(i) e^{-J} |H| Psi >
  !
  ! Dimensions : det_num
  END_DOC

  integer          :: i, j, k, l
  double precision :: T, h_psidet, dij, f, E_noJ, dE

  E_noJ= -0.5d0*psidet_lapl*psidet_inv + E_pot + E_nucl 
  do k=1,det_num
     i = det_coef_matrix_rows(k)
     j = det_coef_matrix_columns(k)
     f = det_alpha_value(i)*det_beta_value(j) * psi_value_inv * jast_value_inv
     ci_dress(k) = E_noJ * f
  enddo

END_PROVIDER

