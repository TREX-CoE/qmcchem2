
! ---

BEGIN_PROVIDER [ double precision, opt_ab_lr, (size_opt_ab_lr) ]

  BEGIN_DOC
  !
  ! matrix elements to optimize alpha and beta of :
  ! Psi = alpha Phi_R e^J + beta Phi_L e^{-J}
  !
  ! Dimensions : 8
  END_DOC

  implicit none
  double precision :: e0, er, el, flr, jj

  e0  = E_pot + E_nucl
  flr = psidet_left_value / psidet_right_value

  er = E_kin_psidet_right + e0 + deltaE_Jmu1b
  el = E_kin_psidet_left  + e0 + deltaE_mJmu1b

  jj = jast_Mu_1b_value * jast_Mu_1b_value

  ! ---

  ! < Phi_R e^{J} | H | Phi_R e^{J} >
  opt_ab_lr(1) = er * jj

  ! < Phi_L e^{-J} | H | Phi_R e^{J} >
  opt_ab_lr(2) = er * flr

  ! < Phi_R e^{J} | H | Phi_L e^{-J} >
  opt_ab_lr(3) = el * flr

  ! < Phi_L e^{-J} | H | Phi_L e^{-J} >
  opt_ab_lr(4) = el * flr * flr / jj

  ! ---

  ! < Phi_R e^{J} | Phi_R e^{J} >
  opt_ab_lr(5) = jj

  ! < Phi_L e^{-J} | Phi_R e^{J} >
  opt_ab_lr(6) = flr

  ! < Phi_R e^{J} | Phi_L e^{-J} >
  opt_ab_lr(7) = flr

  ! < Phi_L e^{-J} | Phi_L e^{-J} >
  opt_ab_lr(8) = flr * flr / jj

  opt_ab_lr_min = min(opt_ab_lr_min, minval(opt_ab_lr))
  opt_ab_lr_max = max(opt_ab_lr_max, maxval(opt_ab_lr))
  SOFT_TOUCH opt_ab_lr_min opt_ab_lr_max
END_PROVIDER

! ---

