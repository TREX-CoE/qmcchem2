
! ---

BEGIN_PROVIDER [ double precision, oo_psik_psi, (size_oo_psik_psi) ]

  BEGIN_DOC
  !
  ! oo_psik_psi(k) = < Psi_k | Psi > / < Psi | Psi >
  !                = < Psi_k / Psi >_{Psi^2}
  ! with:
  !       Psi_k = d Psi / d mo_coef(i,j), k = (i,j)
  !
  ! Dimensions : ao_num*n_oo
  END_DOC

  implicit none
  integer          :: i, j, ii
  double precision :: eps_p, eps_m, eps_f
  double precision :: psi_1, psi_2

  eps_p = 1d-6
  eps_m = 2d-6
  eps_f = 5d+5 * psi_value_inv

  ii = 0
  do j = 1, n_oo
    do i = 1, ao_num

      mo_coef(i,j) = mo_coef(i,j) + eps_p
      TOUCH mo_coef
      psi_1 = psi_value

      mo_coef(i,j) = mo_coef(i,j) - eps_m
      TOUCH mo_coef
      psi_2 = psi_value

      mo_coef(i,j) = mo_coef(i,j) + eps_p
      TOUCH mo_coef

      oo_psik_psi(ii+i) = eps_f * (psi_1 - psi_2) 
    enddo

    ii = ii + ao_num
  enddo 

  oo_psik_psi_min = min(oo_psik_psi_min, minval(oo_psik_psi))
  oo_psik_psi_max = max(oo_psik_psi_max, maxval(oo_psik_psi))
  SOFT_TOUCH oo_psik_psi_min oo_psik_psi_max

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, oo_psik_psil, (size_oo_psik_psil) ]

  BEGIN_DOC
  !
  ! oo_psik_psil(k,l) = < Psi_k | Psi_l > / < Psi | Psi >
  !                   = < [Psi_k / Psi] x [Psi_l / Psi]  >_{Psi^2}
  ! with:
  !       Psi_k = d Psi / d mo_coef(i,j), k = (i,j)
  !
  ! Dimensions : ao_num*n_oo*ao_num*n_oo
  END_DOC

  implicit none
  integer          :: i1, j1, i2, j2, ii
  double precision :: eps_p, eps_m, eps_f, tmp
  double precision :: psi_1, psi_2

  eps_p = 1d-6
  eps_m = 2d-6
  eps_f = 5d+5 * psi_value_inv

  ii = 0
  do j1 = 1, n_oo
    do i1 = 1, ao_num

      mo_coef(i1,j1) = mo_coef(i1,j1) + eps_p
      TOUCH mo_coef
      psi_1 = psi_value

      mo_coef(i1,j1) = mo_coef(i1,j1) - eps_m
      TOUCH mo_coef
      psi_2 = psi_value

      mo_coef(i1,j1) = mo_coef(i1,j1) + eps_p
      TOUCH mo_coef

      tmp = eps_f * (psi_1 - psi_2) 

      do j2 = 1, n_oo
        do i2 = 1, ao_num

          mo_coef(i2,j2) = mo_coef(i2,j2) + eps_p
          TOUCH mo_coef
          psi_1 = psi_value

          mo_coef(i2,j2) = mo_coef(i2,j2) - eps_m
          TOUCH mo_coef
          psi_2 = psi_value

          mo_coef(i2,j2) = mo_coef(i2,j2) + eps_p
          TOUCH mo_coef

          oo_psik_psil(ii+i2) = eps_f * (psi_1 - psi_2) * tmp
        enddo

        ii = ii + ao_num
      enddo
    enddo
  enddo 

  oo_psik_psil_min = min(oo_psik_psil_min, minval(oo_psik_psil))
  oo_psik_psil_max = max(oo_psik_psil_max, maxval(oo_psik_psil))
  SOFT_TOUCH oo_psik_psil_min oo_psik_psil_max

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, oo_psik_H_psi, (size_oo_psik_H_psi) ]

  BEGIN_DOC
  !
  ! oo_psik_H_psi(k) = < Psi_k | H | Psi > / < Psi | Psi >
  !                  = < Psi_k x E_loc / Psi >_{Psi^2}
  ! with:
  !       Psi_k = d Psi / d mo_coef(i,j), k = (i,j)
  !
  ! Dimensions : ao_num*n_oo
  END_DOC

  implicit none
  integer          :: i, j, ii
  double precision :: eps_p, eps_m, eps_f
  double precision :: psi_1, psi_2

  eps_p = 1d-6
  eps_m = 2d-6
  eps_f = 5d+5 * psi_value_inv * E_loc

  ii = 0
  do j = 1, n_oo
    do i = 1, ao_num

      mo_coef(i,j) = mo_coef(i,j) + eps_p
      TOUCH mo_coef
      psi_1 = psi_value

      mo_coef(i,j) = mo_coef(i,j) - eps_m
      TOUCH mo_coef
      psi_2 = psi_value

      mo_coef(i,j) = mo_coef(i,j) + eps_p
      TOUCH mo_coef

      oo_psik_H_psi(ii+i) = eps_f * (psi_1 - psi_2) 
    enddo

    ii = ii + ao_num
  enddo 

  oo_psik_H_psi_min = min(oo_psik_H_psi_min, minval(oo_psik_H_psi))
  oo_psik_H_psi_max = max(oo_psik_H_psi_max, maxval(oo_psik_H_psi))
  SOFT_TOUCH oo_psik_H_psi_min oo_psik_H_psi_max

END_PROVIDER

! ---

