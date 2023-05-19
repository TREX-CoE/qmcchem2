
! ---

 BEGIN_PROVIDER [ logical, use_lr ]
&BEGIN_PROVIDER [ real,    coef_psi_right]
&BEGIN_PROVIDER [ real,    coef_psi_left ]

  BEGIN_DOC
  ! if true, Psi = alpha Phi_R e^J + beta Phi_L e^-J
  END_DOC

  implicit none

  use_lr = .false.
  call get_bi_ortho_mos_use_lr(use_lr)

  if(use_lr) then
    call get_bi_ortho_mos_coef_psi_right(coef_psi_right)
    call get_bi_ortho_mos_coef_psi_left (coef_psi_left )
  else
    coef_psi_right = 0.
    coef_psi_left  = 0.
  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [ logical, inv_sgn_jast ]
&BEGIN_PROVIDER [ double precision, sgn_jast ]

  BEGIN_DOC
  ! +1 for e^{+J}
  ! -1 for e^{-J}
  END_DOC

  implicit none

  sgn_jast = +1.d0

  inv_sgn_jast = .false.
  call get_jastrow_inv_sgn_jast(inv_sgn_jast)
  if(inv_sgn_jast) then
    sgn_jast = -1.d0
  endif

END_PROVIDER

! ---

