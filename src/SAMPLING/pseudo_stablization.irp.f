BEGIN_PROVIDER [ double precision, pseudo_stabilization ]
 implicit none
 BEGIN_DOC
 ! When (W_nl Psi)/Psi < -|(T Psi)/Psi|, we detect a possible 
 ! explosion of the DMC weights. The local energy used in the weights is shifted by 
 ! 1/2 |(W_nl Psi)/Psi| 
 END_DOC
    double precision :: s1, s2
    integer :: i
    real, save :: smin, smax
    logical :: detected

    s1 = 0.d0
    s2 = 0.d0
    detected = .False.
    do i=1,elec_num
      s1 = s1 + e_kin_elec(i)
      s2 = s2 + pseudo_right_non_local(i)
    enddo
    detected = (s2 + dabs(s1) < 0.d0)

    if ( detected ) then
      pseudo_stabilization = -s2 * 0.5d0
    else
      pseudo_stabilization = 0.d0
    endif

END_PROVIDER

