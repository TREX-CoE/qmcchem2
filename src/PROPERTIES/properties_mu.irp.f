
BEGIN_PROVIDER [ double precision, emudiff ]
  implicit none
  BEGIN_DOC
  ! E mu
  END_DOC

  !emudiff = e_loc - energy_mu * jast_value_inv * jast_value_inv
  emudiff = ( e_loc - energy_mu ) * jast_value_inv * jast_value_inv

  emudiff_min = min(emudiff_min,emudiff)
  emudiff_max = max(emudiff_max,emudiff)
  SOFT_TOUCH emudiff_min emudiff_max
END_PROVIDER

BEGIN_PROVIDER [ double precision, Energy_mu ]

  BEGIN_DOC
  ! E mu = < H_mu \Phi / \Phi >_{\Phi^2}
  END_DOC

  implicit none
  integer :: i

  double precision :: lapl
  lapl = 0.d0
  do i=1,elec_num
    lapl += psidet_grad_lapl(4,i)*psidet_inv + jast_elec_mu_lapl(i) + &
    2.d0*psidet_inv * (&
        psidet_grad_lapl(1,i)*jast_elec_mu_grad_x(i) +                   &
        psidet_grad_lapl(2,i)*jast_elec_mu_grad_y(i) +                   &
        psidet_grad_lapl(3,i)*jast_elec_mu_grad_z(i)  ) + ( &
        jast_elec_mu_grad_x(i)*jast_elec_mu_grad_x(i) + &
        jast_elec_mu_grad_y(i)*jast_elec_mu_grad_y(i) + &
        jast_elec_mu_grad_z(i)*jast_elec_mu_grad_z(i) )

  enddo
  Energy_mu = -0.5d0 * lapl + E_nucl + E_pot

  energy_mu_min = min(energy_mu_min,energy_mu)
  energy_mu_max = max(energy_mu_max,energy_mu)
  SOFT_TOUCH energy_mu_min energy_mu_max
END_PROVIDER

BEGIN_PROVIDER [double precision, E_nucl_elec]
 implicit none
 integer :: i,j
 E_nucl_elec = 0.d0
 do i = 1, elec_num
!  E_nucl_elec += E_pot_elec_one(i) + E_pot_elec_two(i)
  E_nucl_elec += E_pot_elec_one(i) 
 enddo
  E_nucl_elec_min = min(E_nucl_elec_min,E_nucl_elec)
  E_nucl_elec_max = max(E_nucl_elec_max,E_nucl_elec)
END_PROVIDER 


 BEGIN_PROVIDER [double precision, Eff_pot_mu_elec, (elec_num)]
&BEGIN_PROVIDER [double precision, Eff_pot_mu_elec_simple, (elec_num)]

  include '../constants.F'

  implicit none
  integer :: i,j
  double precision :: rij, mu

  mu = jast_mu_erf
  Eff_pot_mu_elec = 0.d0


  ! 2body-Jastrow:
  ! 
  ! \Delta_i u_ij + \Delta_j u_ij     = 2 [ (1-erf(mu r_ij))/r_ij - mu exp(-(mu r_ij)^2)/sqrt(pi) ]
  !
  ! (grad_i u_ij)^2 + (grad_j u_ij)^2 = (1-erf(mu r_ij))^2 / 2 

  do i = 1, elec_num
    !DIR$ VECTOR ALIGNED
    !DIR$ LOOP COUNT(50)
    do j = 1, elec_num
      rij = elec_dist(j,i)
      if(i==j)cycle

      Eff_pot_mu_elec(i) = Eff_pot_mu_elec(i) + 0.5d0 * derf(mu * rij) * elec_dist_inv(j,i)
      Eff_pot_mu_elec(i) = Eff_pot_mu_elec(i) + 0.5d0 * mu/dsqpi * dexp(-mu*mu*rij*rij)
      Eff_pot_mu_elec_simple(i) = Eff_pot_mu_elec(i)

      Eff_pot_mu_elec(i) = Eff_pot_mu_elec(i) + 0.5d0 * (-  0.25d0 * (1.d0 - derf(mu*rij))**2.d0 )
                                             
    enddo
  enddo

  ! 1-body Jastrow
  if( jast_1b_type .gt. 0 ) then
    do i = 1, elec_num
      Eff_pot_mu_elec(i) -= 0.5d0 * jast_1b_lapl(i) 
      Eff_pot_mu_elec(i) -= 0.5d0 * jast_1b_grad_sq(i) 
      do j = 1, elec_num
        if(i==j) cycle
        ! + sign for i <--> j
        ! 0.5d0  for double counting
        Eff_pot_mu_elec(i) += 0.5d0 *                                                      &
                            ( grad_j_mu_x(j,i) * ( jast_1b_grad_x(i) - jast_1b_grad_x(j) ) &
                            + grad_j_mu_y(j,i) * ( jast_1b_grad_y(i) - jast_1b_grad_y(j) ) &
                            + grad_j_mu_z(j,i) * ( jast_1b_grad_z(i) - jast_1b_grad_z(j) ) )
      enddo
    enddo
  endif

END_PROVIDER 

 BEGIN_PROVIDER [double precision, Eff_pot_mu ]
 implicit none
 include '../constants.F'
 integer :: i
 Eff_pot_mu = 0.d0
 do i=1,elec_num
  Eff_pot_mu += Eff_pot_mu_elec(i) 
 enddo
  Eff_pot_mu_min = min(Eff_pot_mu_min,Eff_pot_mu)
  Eff_pot_mu_max = max(Eff_pot_mu_max,Eff_pot_mu)
  SOFT_TOUCH Eff_pot_mu_min Eff_pot_mu_max

END_PROVIDER 

BEGIN_PROVIDER [double precision, Eff_pot_mu_simple ]
 implicit none
 include '../constants.F'
 integer :: i
 Eff_pot_mu_simple = 0.d0
 do i=1,elec_num
  Eff_pot_mu_simple += Eff_pot_mu_elec_simple(i) 
 enddo
  Eff_pot_mu_simple_min = min(Eff_pot_mu_simple_min,Eff_pot_mu_simple)
  Eff_pot_mu_simple_max = max(Eff_pot_mu_simple_max,Eff_pot_mu_simple)
  SOFT_TOUCH Eff_pot_mu_simple_min Eff_pot_mu_simple_max

END_PROVIDER 



BEGIN_PROVIDER [double precision, Eff_pot_deriv_mu_elec, (elec_num) ]

  BEGIN_DOC
  !
  ! non-Hermitian term:
  !   - grad_i(tau) . grad_i(\Phi) / \Phi
  !
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: rij, mu

  mu = jast_mu_erf
  Eff_pot_deriv_mu_elec = 0.d0

  ! 2body-Jastrow: (eq A4)
  !   - [ grad_i(tau_mu) . grad_i(\Phi) + grad_j(tau_mu) . grad_j(\Phi) ] / \Phi = 
  !   ( erf(mu r_ij) - 1 ) / ( 2 r_ij \Phi) * [ 
  !                           ( x_i - x_j ) * ( \partial_{x_i} - \partial_{x_j} ) + 
  !                           ( y_i - y_j ) * ( \partial_{y_i} - \partial_{y_j} ) +
  !                           ( z_i - z_j ) * ( \partial_{z_i} - \partial_{z_j} ) ]
  !
  do i = 1, elec_num
    do j = 1, elec_num
      if(i==j)cycle
      rij = elec_dist(i,j)
      Eff_pot_deriv_mu_elec(i) += 0.5d0 * ( derf(mu * rij) - 1.d0 ) * elec_dist_inv(j,i)   &
                                        * ( - elec_dist_vec_x(j,i) * psidet_grad_lapl(1,i) &
                                            - elec_dist_vec_y(j,i) * psidet_grad_lapl(2,i) &
                                            - elec_dist_vec_z(j,i) * psidet_grad_lapl(3,i) ) * psidet_inv
    enddo
  enddo


  ! 1-body Jastrow
  if( jast_1b_type .gt. 0 ) then
    do i = 1, elec_num
      Eff_pot_deriv_mu_elec(i) -= ( jast_1b_grad_x(i) * psidet_grad_lapl(1,i) &
                                  + jast_1b_grad_y(i) * psidet_grad_lapl(2,i) &
                                  + jast_1b_grad_z(i) * psidet_grad_lapl(3,i) ) * psidet_inv
    enddo
  endif

END_PROVIDER 


BEGIN_PROVIDER [double precision, three_body_mu ]
 implicit none
 integer :: i,j,k
 three_body_mu = 0.d0
 do i = 1, elec_num
  do j = i+1, elec_num
   do k = j+1, elec_num
    three_body_mu += grad_j_mu_x(i,j) * grad_j_mu_x(i,k) 
    three_body_mu += grad_j_mu_y(i,j) * grad_j_mu_y(i,k) 
    three_body_mu += grad_j_mu_z(i,j) * grad_j_mu_z(i,k) 

    three_body_mu += grad_j_mu_x(j,i) * grad_j_mu_x(j,k) 
    three_body_mu += grad_j_mu_y(j,i) * grad_j_mu_y(j,k) 
    three_body_mu += grad_j_mu_z(j,i) * grad_j_mu_z(j,k) 

    three_body_mu += grad_j_mu_x(k,i) * grad_j_mu_x(k,j) 
    three_body_mu += grad_j_mu_y(k,i) * grad_j_mu_y(k,j) 
    three_body_mu += grad_j_mu_z(k,i) * grad_j_mu_z(k,j) 
   enddo
  enddo
 enddo
 three_body_mu_min = min(three_body_mu_min,three_body_mu)
 three_body_mu_max = max(three_body_mu_max,three_body_mu)
 SOFT_TOUCH  three_body_mu_min three_body_mu_max
END_PROVIDER 

BEGIN_PROVIDER [double precision, Eff_pot_deriv_mu]
 implicit none
 integer :: i
 Eff_pot_deriv_mu = 0.d0
 do i = 1, elec_num
  Eff_pot_deriv_mu += Eff_pot_deriv_mu_elec(i)
 enddo
  eff_pot_deriv_mu_min = min(eff_pot_deriv_mu_min,eff_pot_deriv_mu)
  eff_pot_deriv_mu_max = max(eff_pot_deriv_mu_max,eff_pot_deriv_mu)
  SOFT_TOUCH eff_pot_deriv_mu_min eff_pot_deriv_mu_max

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ci_dress_mu, (size_ci_dress_mu) ]
  BEGIN_DOC
  ! Dimensions : det_num
  END_DOC
  implicit none
  integer          :: i, j, k, l
  double precision :: T, dij, f, E_noJ, dE
  ! energy_mu = H_mu \Phi / \Phi
  dE = (E_loc - energy_mu) * psi_value_inv * jast_value_inv
  do k = 1, det_num
     i = det_coef_matrix_rows(   k)
     j = det_coef_matrix_columns(k)
     f = det_alpha_value(i) * det_beta_value(j)
     ci_dress_mu(k) = dE * f
  enddo
  ci_dress_mu_min = min(ci_dress_mu_min, minval(ci_dress_mu))
  ci_dress_mu_max = max(ci_dress_mu_max, maxval(ci_dress_mu))
  SOFT_TOUCH ci_dress_mu_min ci_dress_mu_max
END_PROVIDER


BEGIN_PROVIDER [ double precision, ci_dress_mu_opt ]
  BEGIN_DOC
  ! Use for optimizing mu
  END_DOC
  implicit none
  integer          :: i, j, k, l
  double precision :: T, dij, f, E_noJ, dE
  ! energy_mu = H_mu \Phi / \Phi
  dE = (E_loc - energy_mu) * psi_value_inv * jast_value_inv ! PsiJ.J
  k = 1
  i = det_coef_matrix_rows(   k)
  j = det_coef_matrix_columns(k)
  f = det_alpha_value(i) * det_beta_value(j)
  ci_dress_mu_opt = dE * f
  ci_dress_mu_opt = E_loc - energy_mu
  ci_dress_mu_opt_min = min(ci_dress_mu_opt_min, ci_dress_mu_opt)
  ci_dress_mu_opt_max = max(ci_dress_mu_opt_max, ci_dress_mu_opt)
  SOFT_TOUCH ci_dress_mu_opt_min ci_dress_mu_opt_max
END_PROVIDER

