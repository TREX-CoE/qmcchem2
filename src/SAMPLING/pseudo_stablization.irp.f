BEGIN_PROVIDER [ double precision, pseudo_stabilization ]
 implicit none
 BEGIN_DOC
 !
 END_DOC
    double precision :: s1, s2
    integer :: i

    s1 = 0.d0
    s2 = 0.d0
    do i=1,elec_num
      s1 = s1 + v_pseudo_local(i) + pseudo_right_non_local(i)
      s2 = s2 + pseudo_right_non_local(i)
    enddo

    if ( (s1<0.d0).and.(s2<0.d0) ) then
      pseudo_stabilization = -s2
    else
      pseudo_stabilization = 0.d0
    endif

END_PROVIDER


