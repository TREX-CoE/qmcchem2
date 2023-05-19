
! ---

BEGIN_PROVIDER [ double precision, H_DJ, (det_num) ]

  implicit none
  integer          :: l, m, n, e, ee
  double precision :: g, h, T, V, j_lapl_inv

  ! (Lapl J)/J
  j_lapl_inv = jast_lapl_jast_inv_tot

  do l = 1, det_num
    m  = det_coef_matrix_rows   (l)
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
        V += det_right_alpha_value(m) * det_right_beta_pseudo(e,n)
      enddo
    endif

    H_DJ(l) = ( -0.5d0*T + V ) * psidet_right_inv
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, jast_lapl_jast_inv_tot]

  BEGIN_DOC
  ! (Lapl J)/J
  END_DOC

  implicit none
  integer          :: e
  double precision :: j_tmp

  j_tmp = 0.d0
  do e = 1, elec_num
    j_tmp = j_tmp + jast_lapl_jast_inv(e)
  enddo
  jast_lapl_jast_inv_tot = j_tmp

END_PROVIDER

! ---

