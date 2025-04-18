BEGIN_PROVIDER [ double precision, pseudo_stabilization ]
 implicit none
 BEGIN_DOC
 ! Eq(33) of 10.1063/1.465195
 ! When the single-electron drift differs too much from
 ! (-1+sqrt(1+2V^2 tau))/(V^2 tau)
 ! The DMC weight is set to 1.
 END_DOC
    double precision :: s1, s2
    integer :: i
    logical :: detected

    s1 = 0.d0
    s2 = 0.d0
    do i=1,elec_num
      s1 = s1 + e_kin_elec(i) * drift_renorm(i)
      s2 = s2 + pseudo_right_non_local(i)
    enddo
    detected = ( (s2 < 0.d0) .and. dabs(s2) > dabs(s1) )

    if ( detected ) then
      pseudo_stabilization = e_ref - e_loc
    else
      pseudo_stabilization = 0.d0
    endif

END_PROVIDER

