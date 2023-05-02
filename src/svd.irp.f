
! ---

 BEGIN_PROVIDER [ integer, n_svd_alpha ]
&BEGIN_PROVIDER [ integer, n_svd_beta ]
&BEGIN_PROVIDER [ integer, n_svd_coefs  ]
&BEGIN_PROVIDER [ logical, use_svd ]

  BEGIN_DOC
  ! U: n_det_alpha_num x n_svd_alpha
  ! V: n_det_beta_num  x n_svd_beta
  END_DOC

  implicit none

  n_svd_alpha = det_alpha_num
  n_svd_beta  = det_beta_num
  call get_spindeterminants_n_svd_alpha(n_svd_alpha)
  call get_spindeterminants_n_svd_beta (n_svd_beta )

  n_svd_coefs = min(n_svd_alpha, n_svd_beta)

  use_svd = .false.

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, psi_svd_alpha, (det_alpha_num, n_svd_alpha, n_states) ]
&BEGIN_PROVIDER [ double precision, psi_svd_beta , (det_beta_num , n_svd_beta , n_states) ]

   BEGIN_DOC
   ! U = psi_svd_alpha: det_alpha_num x n_svd_alpha
   ! V = psi_det_beta : det_beta_num  x n_svd_beta
   END_DOC

   implicit none
   integer :: i

   !call get_spindeterminants_psi_svd_alpha(psi_svd_alpha)
   !call get_spindeterminants_psi_svd_beta (psi_svd_beta )

   open( unit = 11, form = "unformatted" &
       , file = "/users/p22001/p22001aa/tmpdir/dress_opt/QMC/cr2_work/Jmu/IT1/SVD_basis/P3/U_unf" &
       , action = "read" )
     do i = 1, n_svd_alpha
       read(11) psi_svd_alpha(:,i,1)
     enddo
   close(11)

   open( unit = 11, form = "unformatted" &
       , file = "/users/p22001/p22001aa/tmpdir/dress_opt/QMC/cr2_work/Jmu/IT1/SVD_basis/P3/V_unf" &
       , action = "read" )
     do i = 1, n_svd_beta
       read(11) psi_svd_beta(:,i,1)
     enddo
   close(11)


END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, psi_svd_coefs, (n_svd_coefs, n_states) ]

   implicit none

   call get_spindeterminants_psi_svd_coefs(psi_svd_coefs)

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, det_alpha_value_svd, (n_svd_alpha) ]
&BEGIN_PROVIDER [ double precision, det_beta_value_svd , (n_svd_beta ) ]

  implicit none

  ! det_alpha_value_svd = psi_svd_alpha.T @ det_right_alpha_value
  call dgemv('T', det_alpha_num, n_svd_alpha                                &
    , 1.d0, psi_svd_alpha, size(psi_svd_alpha, 1), det_right_alpha_value, 1 &
    , 0.d0, det_alpha_value_svd, 1)

  ! det_beta_value_svd = psi_svd_beta.T @ det_right_beta_value
  call dgemv('T', det_beta_num, n_svd_beta                               &
    , 1.d0, psi_svd_beta, size(psi_svd_beta, 1), det_right_beta_value, 1 &
    , 0.d0, det_beta_value_svd, 1)

END_PROVIDER

! ---

 BEGIN_PROVIDER [ integer, la   ]
&BEGIN_PROVIDER [ integer, lb   ]
&BEGIN_PROVIDER [ integer, ld   ]
&BEGIN_PROVIDER [ integer, lla  ]
&BEGIN_PROVIDER [ integer, llb  ]
&BEGIN_PROVIDER [ integer, lld  ]
&BEGIN_PROVIDER [ integer, n_dress_svd ]

   BEGIN_DOC
   ! svd 4b dim for dressing
   ! lx : min along x, with x = a (row), b (col), d (diag)
   ! llx: max along x, with x = a (row), b (col), d (diag)
   END_DOC

   implicit none

   call get_dmc_dress_la(la)
   call get_dmc_dress_lb(lb)
   call get_dmc_dress_ld(ld)
   call get_dmc_dress_lla(lla)
   call get_dmc_dress_llb(llb)
   call get_dmc_dress_lld(lld)

   ! dim of delta_svd vector
   n_dress_svd = ( ld * ld       &
                 + lld - ld      &
                 + (lla-ld) * lb &
                 + (llb-ld) * la )

END_PROVIDER

! ---




