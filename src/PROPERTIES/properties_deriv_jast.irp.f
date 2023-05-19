
BEGIN_PROVIDER  [ double precision, E_deriv_nucPar_loc1, (size_E_deriv_nucPar_loc1) ]
  implicit none
  BEGIN_DOC
  ! Local energy variation  with respect to nuclear parameters
  ! 
  ! Dimensions : nucl_num
  END_DOC

  integer :: i

  do i = 1, nucl_num
      E_deriv_nucPar_loc1(i) = E_loc*jast_elec_Simple_deriv_nucPar(i)
  enddo

  E_deriv_nucPar_loc1_min = min(E_deriv_nucPar_loc1_min,minval(E_deriv_nucPar_loc1))
  E_deriv_nucPar_loc1_max = max(E_deriv_nucPar_loc1_max,maxval(E_deriv_nucPar_loc1))
  SOFT_TOUCH E_deriv_nucPar_loc1_min E_deriv_nucPar_loc1_max

END_PROVIDER


BEGIN_PROVIDER  [ double precision, E_deriv_nucPar_loc2, (size_E_deriv_nucPar_loc2) ]
  implicit none
  BEGIN_DOC
  ! Local energy variation  with respect to nuclear parameters
  ! 
  ! Dimensions : nucl_num
  END_DOC

  integer :: i

  do i=1,nucl_num
      E_deriv_nucPar_loc2(i) = jast_elec_Simple_deriv_nucPar(i)
  enddo

  E_deriv_nucPar_loc2_min = min(E_deriv_nucPar_loc2_min,minval(E_deriv_nucPar_loc2))
  E_deriv_nucPar_loc2_max = max(E_deriv_nucPar_loc2_max,maxval(E_deriv_nucPar_loc2))
  SOFT_TOUCH E_deriv_nucPar_loc2_min E_deriv_nucPar_loc2_max
END_PROVIDER



BEGIN_PROVIDER  [ double precision, E_deriv_bPar_loc1 ]
  implicit none
  BEGIN_DOC
  ! Local energy variation  with respect to parameter b
  END_DOC

  E_deriv_bPar_loc1 = E_loc*J_deriv_bPar_ex

  E_deriv_bPar_loc1_min = min(E_deriv_bPar_loc1_min, E_deriv_bPar_loc1)
  E_deriv_bPar_loc1_max = max(E_deriv_bPar_loc1_max, E_deriv_bPar_loc1)
  SOFT_TOUCH E_deriv_bPar_loc1_min E_deriv_bPar_loc1_max
END_PROVIDER

BEGIN_PROVIDER  [ double precision, E_deriv_bPar_loc2 ]
  implicit none
  BEGIN_DOC
  ! Local energy variation  with respect to parameter b
  END_DOC

  E_deriv_bPar_loc2 = J_deriv_bPar_ex

  E_deriv_bPar_loc2_min = min(E_deriv_bPar_loc2_min,  E_deriv_bPar_loc2)
  E_deriv_bPar_loc2_max = max(E_deriv_bPar_loc2_max,  E_deriv_bPar_loc2)
  SOFT_TOUCH E_deriv_bPar_loc2_min E_deriv_bPar_loc2_max
END_PROVIDER





BEGIN_PROVIDER  [ double precision, J_deriv_nucPar_verif, (size_J_deriv_nucPar_verif) ]
  implicit none
  BEGIN_DOC
  ! Jastrow variation  with respect to nuclear parameters
  ! 
  ! Dimensions : nucl_num
  END_DOC

  integer          :: i, j
  double precision :: eps = 1d-7, der1, der2, pos0

  do j = 1, nucl_num
    !!!
    pos0 = eps * jast_pen(j)
    !!!
    jast_pen(j) = jast_pen(j) + pos0
    TOUCH jast_pen
    der1 = 0.d0
    !DIR$ LOOP COUNT (100)
    do i = 1, elec_num
      der1 += jast_elec_Simple_value(i)
    end do
    !!!
    jast_pen(j) = jast_pen(j) - 2.d0 * pos0
    TOUCH jast_pen
    der2 = 0.d0
    !DIR$ LOOP COUNT (100)
    do i = 1, elec_num
      der2 += jast_elec_Simple_value(i)
    end do
    !!!
    J_deriv_nucPar_verif(j) = 0.5d0 * (der1 - der2) / pos0
    !!!
  end do

  J_deriv_nucPar_verif_min = min(J_deriv_nucPar_verif_min,minval(J_deriv_nucPar_verif))
  J_deriv_nucPar_verif_max = max(J_deriv_nucPar_verif_max,maxval(J_deriv_nucPar_verif))
  SOFT_TOUCH J_deriv_nucPar_verif_min J_deriv_nucPar_verif_max

END_PROVIDER



BEGIN_PROVIDER  [ double precision, J_deriv_nucPar_ex, (size_J_deriv_nucPar_ex) ]
  implicit none
  BEGIN_DOC
  ! Jastrow variation  with respect to nuclear parameters
  ! 
  ! Dimensions : nucl_num
  END_DOC

  integer          :: j

  do j = 1, nucl_num
    J_deriv_nucPar_ex(j) = jast_elec_Simple_deriv_nucPar(j)
  end do

  J_deriv_nucPar_ex_min = min(J_deriv_nucPar_ex_min,minval(J_deriv_nucPar_ex))
  J_deriv_nucPar_ex_max = max(J_deriv_nucPar_ex_max,maxval(J_deriv_nucPar_ex))
  SOFT_TOUCH J_deriv_nucPar_ex_min J_deriv_nucPar_ex_max

END_PROVIDER



BEGIN_PROVIDER  [ double precision, J_deriv_bPar_ex ]
  implicit none
  BEGIN_DOC
  ! Jastrow variation with respect to parameter b
  ! we suppose that jast_b_up_up = jast_b_up_dn
  END_DOC

  double precision :: tmp, Jtmp, a, b, rij
  integer          :: i, j

  b = jast_b_up_up ! and also = jast_b_up_dn

  ! parallele spin up-up
  Jtmp = 0.d0
  a    = jast_a_up_up
  do j = 1, elec_alpha_num-1
    !DIR$ LOOP COUNT (50)
    do i = j+1, elec_alpha_num
      rij  = elec_dist(i,j)
      tmp  = rij/(1.d0+b*rij) 
      Jtmp = Jtmp + tmp * tmp
    enddo
  enddo
  do j = elec_alpha_num+1, elec_num-1
    !DIR$ LOOP COUNT (50)
    do i = j+1, elec_num
      rij  = elec_dist(i,j)
      tmp  = rij/(1.d0+b*rij) 
      Jtmp = Jtmp + tmp * tmp
    enddo
  enddo
  J_deriv_bPar_ex = -1.d0 * a * Jtmp



  ! anti-parallele spin
  Jtmp = 0.d0
  a    = jast_a_up_dn
  do j = 1, elec_alpha_num
  !DIR$ LOOP COUNT (50)
    do i = elec_alpha_num+1, elec_num
      rij  = elec_dist(i,j)
      tmp  = rij/(1.d0+b*rij) 
      Jtmp = Jtmp + tmp * tmp
    enddo
  enddo
  J_deriv_bPar_ex = J_deriv_bPar_ex - a * Jtmp


  J_deriv_bPar_ex_min = min(J_deriv_bPar_ex_min, J_deriv_bPar_ex)
  J_deriv_bPar_ex_max = max(J_deriv_bPar_ex_max, J_deriv_bPar_ex)
  SOFT_TOUCH J_deriv_bPar_ex_min J_deriv_bPar_ex_max

END_PROVIDER



BEGIN_PROVIDER  [ double precision, J_deriv_bPar_verif ]
  implicit none
  BEGIN_DOC
  ! verification: Jastrow variation with respect to parameter b
  ! we suppose that jast_b_up_up = jast_b_up_dn
  END_DOC

  integer          :: i
  double precision :: eps = 1d-7, der1, der2, pos0

  pos0 = eps * jast_b_up_up
  ! !!!
  jast_b_up_up = jast_b_up_up + pos0
  TOUCH jast_b_up_up
  jast_b_up_dn = jast_b_up_up
  TOUCH jast_b_up_dn
  der1 = 0.d0
  !DIR$ LOOP COUNT (100)
  do i = 1, elec_num
    der1 += jast_elec_Simple_value(i)
  end do
  ! !!!
  jast_b_up_up = jast_b_up_up - 2.d0 * pos0
  TOUCH jast_b_up_up
  jast_b_up_dn = jast_b_up_up
  TOUCH jast_b_up_dn
  der2 = 0.d0
  !DIR$ LOOP COUNT (100)
  do i = 1, elec_num
    der2 += jast_elec_Simple_value(i)
  end do
  ! !!!
  J_deriv_bPar_verif = 0.5d0 * (der1 - der2) / pos0


  J_deriv_bPar_verif_min = min(J_deriv_bPar_verif_min, J_deriv_bPar_verif)
  J_deriv_bPar_verif_max = max(J_deriv_bPar_verif_max, J_deriv_bPar_verif)
  SOFT_TOUCH J_deriv_bPar_verif_min J_deriv_bPar_verif_max

END_PROVIDER
