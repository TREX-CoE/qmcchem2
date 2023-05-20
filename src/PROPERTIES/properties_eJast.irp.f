! ---

! -----------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, Eloc_noJ ]

  BEGIN_DOC
  ! Eloc_noJ = H Phi / Phi 
  END_DOC

  implicit none

  !Eloc_noJ = -0.5d0 * psidet_right_lapl * psidet_right_inv + E_pot + E_nucl
  Eloc_noJ = E_kin_psidet_right + E_pot + E_nucl

  Eloc_noJ_min = min(Eloc_noJ_min, Eloc_noJ)
  Eloc_noJ_max = max(Eloc_noJ_max, Eloc_noJ)
  SOFT_TOUCH Eloc_noJ_min Eloc_noJ_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_noJ_erf ]

  BEGIN_DOC
  ! Eloc_noJ_erf = H_erf Phi / Phi 
  END_DOC

  implicit none

  Eloc_noJ_erf = E_kin_psidet_right + E_pot_erf + E_nucl

  Eloc_noJ_erf_min = min(Eloc_noJ_erf_min, Eloc_noJ_erf)
  Eloc_noJ_erf_max = max(Eloc_noJ_erf_max, Eloc_noJ_erf)
  SOFT_TOUCH Eloc_noJ_erf_min Eloc_noJ_erf_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_Jsimple ]

  BEGIN_DOC
  ! Eloc_Jsimple = H ( Phi * e^{Jsimple} ) / ( Phi * e^{Jsimple} ) 
  END_DOC

  implicit none

  Eloc_Jsimple = Eloc_noJ + deltaE_Jsimple 

  Eloc_Jsimple_min = min(Eloc_Jsimple_min, Eloc_Jsimple)
  Eloc_Jsimple_max = max(Eloc_Jsimple_max, Eloc_Jsimple)
  SOFT_TOUCH Eloc_Jsimple_min Eloc_Jsimple_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_Jmu1b ]

  BEGIN_DOC
  ! Eloc_Jmu1b = H ( Phi * e^{Jmu1b} ) / ( Phi * e^{Jmu1b} ) 
  END_DOC

  implicit none

  Eloc_Jmu1b = Eloc_noJ + deltaE_Jmu1b

  Eloc_Jmu1b_min = min(Eloc_Jmu1b_min, Eloc_Jmu1b)
  Eloc_Jmu1b_max = max(Eloc_Jmu1b_max, Eloc_Jmu1b)
  SOFT_TOUCH Eloc_Jmu1b_min Eloc_Jmu1b_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_Jmu1b_no3b ]

  BEGIN_DOC
  ! Eloc_Jmu1b_no3b = H ( Phi * e^{Jmu1b} ) / ( Phi * e^{Jmu1b} ) 
  !                 - deltaE_Jmu_grad_3b 
  END_DOC

  implicit none

  Eloc_Jmu1b_no3b = Eloc_noJ + deltaE_Jmu1b_no3b

  Eloc_Jmu1b_no3b_min = min(Eloc_Jmu1b_no3b_min, Eloc_Jmu1b_no3b)
  Eloc_Jmu1b_no3b_max = max(Eloc_Jmu1b_no3b_max, Eloc_Jmu1b_no3b)
  SOFT_TOUCH Eloc_Jmu1b_no3b_min Eloc_Jmu1b_no3b_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_Jpsi ]

  BEGIN_DOC
  ! Eloc_Jpsi = H ( Phi * e^{Jpsi} ) / ( Phi * e^{Jpsi} ) 
  END_DOC

  implicit none

  Eloc_Jpsi = Eloc_noJ + deltaE_Jpsi

  Eloc_Jpsi_min = min(Eloc_Jpsi_min, Eloc_Jpsi)
  Eloc_Jpsi_max = max(Eloc_Jpsi_max, Eloc_Jpsi)
  SOFT_TOUCH Eloc_Jpsi_min Eloc_Jpsi_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_mJpsi ]

  BEGIN_DOC
  ! Eloc_mJpsi = H ( Phi * e^{-Jpsi} ) / ( Phi * e^{-Jpsi} ) 
  END_DOC

  implicit none

  Eloc_mJpsi = Eloc_noJ + deltaE_mJpsi

  Eloc_mJpsi_min = min(Eloc_mJpsi_min, Eloc_mJpsi)
  Eloc_mJpsi_max = max(Eloc_mJpsi_max, Eloc_mJpsi)
  SOFT_TOUCH Eloc_mJpsi_min Eloc_mJpsi_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_Jmu1b_lr ]

  BEGIN_DOC
  ! Eloc_Jmu1b_lr = ( Phi_left / Phi_right ) x E_right 
  ! with :
  ! E_right = H ( Phi_right * e^{Jmu1b} ) / ( Phi_right * e^{Jmu1b} ) 
  END_DOC

  implicit none

  Eloc_Jmu1b_lr = (Eloc_noJ + deltaE_Jmu1b) * psidet_left_value / psidet_right_value

  Eloc_Jmu1b_lr_min = min(Eloc_Jmu1b_lr_min, Eloc_Jmu1b_lr)
  Eloc_Jmu1b_lr_max = max(Eloc_Jmu1b_lr_max, Eloc_Jmu1b_lr)
  SOFT_TOUCH Eloc_Jmu1b_lr_min Eloc_Jmu1b_lr_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Eloc_Jmu1b_lr_no3b ]

  BEGIN_DOC
  ! Eloc_Jmu1b_lr_no3b = ( Phi_left / Phi_right ) x [E_right - E_right_3b]
  END_DOC

  implicit none

  Eloc_Jmu1b_lr_no3b = (Eloc_noJ + deltaE_Jmu1b_no3b) * psidet_left_value / psidet_right_value

  Eloc_Jmu1b_lr_no3b_min = min(Eloc_Jmu1b_lr_no3b_min, Eloc_Jmu1b_lr_no3b)
  Eloc_Jmu1b_lr_no3b_max = max(Eloc_Jmu1b_lr_no3b_max, Eloc_Jmu1b_lr_no3b)
  SOFT_TOUCH Eloc_Jmu1b_lr_no3b_min Eloc_Jmu1b_lr_no3b_max
END_PROVIDER

! -----------------------------------------------------------------------------

! ---

! -----------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, deltaE_Jpsi ]

  BEGIN_DOC
  ! deltaE_Jpsi = H ( Phi * e^{Jpsi} ) / ( Phi * e^{Jpsi} ) 
  !             - H ( Phi ) / Phi 
  END_DOC

  implicit none

  deltaE_Jpsi = deltaE_Jpsi_lapl + deltaE_Jpsi_nonh + deltaE_Jpsi_grad

  deltaE_Jpsi_min = min(deltaE_Jpsi_min, deltaE_Jpsi)
  deltaE_Jpsi_max = max(deltaE_Jpsi_max, deltaE_Jpsi)
  SOFT_TOUCH deltaE_Jpsi_min deltaE_Jpsi_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_mJpsi ]

  BEGIN_DOC
  ! deltaE_mJpsi = H ( Phi * e^{-Jpsi} ) / ( Phi * e^{-Jpsi} ) 
  !              - H ( Phi ) / Phi 
  END_DOC

  implicit none

  deltaE_mJpsi = - deltaE_Jpsi_lapl - deltaE_Jpsi_nonh + deltaE_Jpsi_grad

  deltaE_mJpsi_min = min(deltaE_mJpsi_min, deltaE_mJpsi)
  deltaE_mJpsi_max = max(deltaE_mJpsi_max, deltaE_mJpsi)
  SOFT_TOUCH deltaE_mJpsi_min deltaE_mJpsi_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jsimple ]

  BEGIN_DOC
  ! deltaE_Jsimple = H ( Phi * e^{Jsimple} ) / ( Phi * e^{Jsimple} ) 
  !                - H ( Phi ) / Phi 
  END_DOC

  implicit none

  deltaE_Jsimple = deltaE_Jsimple_lapl + deltaE_Jsimple_nonh + deltaE_Jsimple_grad

  deltaE_Jsimple_min = min(deltaE_Jsimple_min, deltaE_Jsimple)
  deltaE_Jsimple_max = max(deltaE_Jsimple_max, deltaE_Jsimple)
  SOFT_TOUCH deltaE_Jsimple_min deltaE_Jsimple_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jsimple_reg ]

  BEGIN_DOC
  ! deltaE_Jsimple_reg = H ( Phi * e^{Jsimple} ) / ( Phi * e^{Jsimple} ) 
  !                    - H ( Phi ) / Phi + 1 / r
  END_DOC

  implicit none

  deltaE_Jsimple_reg = deltaE_Jsimple_lapl_reg + deltaE_Jsimple_nonh + deltaE_Jsimple_grad

  deltaE_Jsimple_reg_min = min(deltaE_Jsimple_reg_min, deltaE_Jsimple_reg)
  deltaE_Jsimple_reg_max = max(deltaE_Jsimple_reg_max, deltaE_Jsimple_reg)
  SOFT_TOUCH deltaE_Jsimple_reg_min deltaE_Jsimple_reg_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b ]

  BEGIN_DOC
  ! deltaE_Jmu1b = H ( Phi * e^{Jmu1b} ) / ( Phi * e^{Jmu1b} ) 
  !              - H ( Phi ) / Phi 
  END_DOC

  implicit none

  deltaE_Jmu1b = deltaE_Jmu1b_lapl + deltaE_Jmu1b_nonh + deltaE_Jmu1b_grad

  deltaE_Jmu1b_min = min(deltaE_Jmu1b_min, deltaE_Jmu1b)
  deltaE_Jmu1b_max = max(deltaE_Jmu1b_max, deltaE_Jmu1b)
  SOFT_TOUCH deltaE_Jmu1b_min deltaE_Jmu1b_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_reg ]

  BEGIN_DOC
  ! deltaE_Jmu1b_reg = H ( Phi * e^{Jmu1b} ) / ( Phi * e^{Jmu1b} ) 
  !                  - H ( Phi ) / Phi + 1 / r
  END_DOC

  implicit none

  deltaE_Jmu1b_reg = deltaE_Jmu1b_lapl_reg + deltaE_Jmu1b_nonh + deltaE_Jmu1b_grad

  deltaE_Jmu1b_reg_min = min(deltaE_Jmu1b_reg_min, deltaE_Jmu1b_reg)
  deltaE_Jmu1b_reg_max = max(deltaE_Jmu1b_reg_max, deltaE_Jmu1b_reg)
  SOFT_TOUCH deltaE_Jmu1b_reg_min deltaE_Jmu1b_reg_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_no3b ]

  BEGIN_DOC
  ! deltaE_Jmu1b_no3b = H ( Phi * e^{Jmu1b} ) / ( Phi * e^{Jmu1b} ) 
  !                   - deltaE_Jmu_grad_3b
  !                   - H ( Phi ) / Phi 
  END_DOC

  implicit none

  deltaE_Jmu1b_no3b = deltaE_Jmu1b_lapl + deltaE_Jmu1b_nonh + deltaE_Jmu1b_grad_no3b

  deltaE_Jmu1b_no3b_min = min(deltaE_Jmu1b_no3b_min, deltaE_Jmu1b_no3b)
  deltaE_Jmu1b_no3b_max = max(deltaE_Jmu1b_no3b_max, deltaE_Jmu1b_no3b)
  SOFT_TOUCH deltaE_Jmu1b_no3b_min deltaE_Jmu1b_no3b_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_Jmu1b_no3b_reg ]

  BEGIN_DOC
  ! deltaE_Jmu1b_no3b_reg = H ( Phi * e^{Jmu1b} ) / ( Phi * e^{Jmu1b} ) 
  !                       - deltaE_Jmu_grad_3b
  !                       - H ( Phi ) / Phi + 1 / r
  END_DOC

  implicit none

  deltaE_Jmu1b_no3b_reg = deltaE_Jmu1b_lapl_reg + deltaE_Jmu1b_nonh + deltaE_Jmu1b_grad_no3b

  deltaE_Jmu1b_no3b_reg_min = min(deltaE_Jmu1b_no3b_reg_min, deltaE_Jmu1b_no3b_reg)
  deltaE_Jmu1b_no3b_reg_max = max(deltaE_Jmu1b_no3b_reg_max, deltaE_Jmu1b_no3b_reg)
  SOFT_TOUCH deltaE_Jmu1b_no3b_reg_min deltaE_Jmu1b_no3b_reg_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, deltaE_mJmu1b ]

  BEGIN_DOC
  ! deltaE_mJmu1b = H ( Phi * e^{-Jmu1b} ) / ( Phi * e^{-Jmu1b} ) 
  !               - H ( Phi ) / Phi 
  END_DOC

  implicit none

  deltaE_mJmu1b = - deltaE_Jmu1b_lapl - deltaE_mJmu1b_nonh + deltaE_Jmu1b_grad

  deltaE_mJmu1b_min = min(deltaE_mJmu1b_min, deltaE_mJmu1b)
  deltaE_mJmu1b_max = max(deltaE_mJmu1b_max, deltaE_mJmu1b)
  SOFT_TOUCH deltaE_mJmu1b_min deltaE_mJmu1b_max
END_PROVIDER

! ---

! -----------------------------------------------------------------------------

! ---

! -----------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, EJsimple_m_EJmu1b ]

  BEGIN_DOC
  !
  ! EJsimple_m_EJmu1b = Eloc_Jsimple       - Eloc_Jmu1b       
  !                   = deltaE_Jsimple     - deltaE_Jmu1b     
  !                   = deltaE_Jsimple_reg - deltaE_Jmu1b_reg 
  ! with:
  !        Eloc_Jsimple = H (Phi e^{Jsimple}) / (Phi e^{Jsimple})
  !                     = Eloc_noJ + deltaE_Jsimple
  !        Eloc_Jmu1b   = H (Phi e^{Jmu1b}  ) / (Phi e^{Jmu1b}  )
  !                     = Eloc_noJ + deltaE_Jmu1b
  ! and:
  !        Eloc_Jsimple_reg = Eloc_Jsimple + 1 / r
  !        Eloc_Jmu1b_reg   = Eloc_Jmu1b   + 1 / r
  !
  END_DOC

  implicit none

  !EJsimple_m_EJmu1b = deltaE_Jsimple - deltaE_Jmu1b
  EJsimple_m_EJmu1b = deltaE_Jsimple_reg - deltaE_Jmu1b_reg 

  EJsimple_m_EJmu1b_min = min(EJsimple_m_EJmu1b_min, EJsimple_m_EJmu1b)
  EJsimple_m_EJmu1b_max = max(EJsimple_m_EJmu1b_max, EJsimple_m_EJmu1b)
  SOFT_TOUCH EJsimple_m_EJmu1b_min EJsimple_m_EJmu1b_max
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, EJsimple_m_EJmu1b_no3b ]

  BEGIN_DOC
  !  = Eloc_Jsimple       - Eloc_Jmu1b_no3b
  !  = deltaE_Jsimple     - deltaE_Jmu1b_no3b  
  !  = deltaE_Jsimple_reg - deltaE_Jmu1b_no3b_reg
  END_DOC

  implicit none

  !EJsimple_m_EJmu1b_no3b = deltaE_Jsimple - deltaE_Jmu1b_no3b
  EJsimple_m_EJmu1b_no3b = deltaE_Jsimple_reg - deltaE_Jmu1b_no3b_reg

  EJsimple_m_EJmu1b_no3b_min = min(EJsimple_m_EJmu1b_no3b_min, EJsimple_m_EJmu1b_no3b)
  EJsimple_m_EJmu1b_no3b_max = max(EJsimple_m_EJmu1b_no3b_max, EJsimple_m_EJmu1b_no3b)
  SOFT_TOUCH EJsimple_m_EJmu1b_no3b_min EJsimple_m_EJmu1b_no3b_max
END_PROVIDER

! -----------------------------------------------------------------------------

! ---

! -----------------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, EJsimple_m_Eerf ]

  BEGIN_DOC
  !
  ! EJsimple_m_Eerf = Eloc_Jsimple - Eloc_noJ_erf
  !
  ! with:
  !        Eloc_Jsimple = H (Phi e^{Jsimple}) / (Phi e^{Jsimple})
  !        Eloc_noJ_erf = H_erf (Phi)         / (Phi)
  !
  END_DOC

  implicit none

  EJsimple_m_Eerf = Eloc_Jsimple - Eloc_noJ_erf

  EJsimple_m_Eerf_min = min(EJsimple_m_Eerf_min, EJsimple_m_Eerf)
  EJsimple_m_Eerf_max = max(EJsimple_m_Eerf_max, EJsimple_m_Eerf)
  SOFT_TOUCH EJsimple_m_Eerf_min EJsimple_m_Eerf_max
END_PROVIDER

! -----------------------------------------------------------------------------


