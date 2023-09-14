



BEGIN_PROVIDER [ double precision, eHj_m_eH ]

  BEGIN_DOC
  ! < [ E_loc - Eloc_noJ * e^{-2J} ] >_{Psi^2}
  END_DOC

  implicit none

  eHj_m_eH = E_loc - Eloc_noJ * jast_value_inv * jast_value_inv 

  eHj_m_eH_min = min(eHj_m_eH_min, eHj_m_eH)
  eHj_m_eH_max = max(eHj_m_eH_max, eHj_m_eH)
  SOFT_TOUCH eHj_m_eH_min eHj_m_eH_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision , jast_elec_check, (elec_num_8) ]

  implicit none
  integer          :: i, j
  double precision :: rij, tmp

  do i = 1, elec_num
    jast_elec_check(i) = 0.d0
  enddo

  !do j = 1, elec_alpha_num
  !  !DIR$ LOOP COUNT (50)
  !  do i = 1, j-1
  !    rij = elec_dist(i,j)
  !    tmp = 1.d0/(elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
  !    jast_elec_check(i) += tmp
  !    jast_elec_check(j) += tmp
  !  enddo
  !enddo
  !do j = elec_alpha_num+1, elec_num
  !  !DIR$ LOOP COUNT (100)
  !  do i = j+1, elec_num
  !    rij = elec_dist(i,j)
  !    tmp = 1.d0/(elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
  !    jast_elec_check(i) += tmp
  !    jast_elec_check(j) += tmp
  !  enddo
  !enddo
  !do j = 1, elec_alpha_num
  !  !DIR$ LOOP COUNT (100)
  !  do i = elec_alpha_num+1, elec_num
  !    rij = elec_dist(i,j)
  !    tmp = 1.d0 / (elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
  !    jast_elec_check(i) += tmp
  !    jast_elec_check(j) += tmp
  !  enddo
  !enddo

  !do i = 1, elec_num
  !  do j = 1, elec_num
  !    if(i==j) cycle
  !    rij = elec_dist(i,j)
  !    tmp = 1.d0 / (elec_dist(i,j)*(1.d0+rij*(3.d0+rij*(3.d0+rij))))
  !    jast_elec_check(i) -= tmp
  !  enddo
  !enddo

  do i = 1, elec_num
    jast_elec_check(i) = jast_elec_Simple_lapl_reg(i) - jast_elec_Simple_lapl(i)
    do j = 1, elec_num
      if(i==j) cycle
      rij = elec_dist(i,j)
      jast_elec_check(i) += 1.d0/rij
    enddo
  enddo


  jast_elec_check_min = min(jast_elec_check_min, minval(jast_elec_check))
  jast_elec_check_max = max(jast_elec_check_max, maxval(jast_elec_check))
  SOFT_TOUCH jast_elec_check_min jast_elec_check_max 
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_check ]

  implicit none

  deltaE_check = (deltaE_Jmu_grad_3b + deltaE_Jmu1b_grad_no3b) - deltaE_Jmu1b_grad
  !deltaE_check = (EJsimple_m_EJmu1b_no3b - deltaE_Jmu_grad_3b) - EJsimple_m_EJmu1b

  deltaE_check_min = min(deltaE_check_min, deltaE_check)
  deltaE_check_max = max(deltaE_check_max, deltaE_check)
  SOFT_TOUCH deltaE_check_min deltaE_check_max 
END_PROVIDER

! ---
