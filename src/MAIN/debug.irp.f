program debug
 implicit none

 logical :: test_lapl    = .False.
 logical :: test_mo      = .False.
 logical :: test_E_loc   = .False.
 logical :: test_psi     = .False.
 logical :: test_E_kin   = .False.
 logical :: test_OrbJast_value = .True. 
 logical :: test_OrbJast_grad  = .False.  
 logical :: test_OrbJast_lapl  = .False. 
 logical :: test_OrbJast_psi_value = .False.
 logical :: test_OrbJast_psi_grad  = .False.
 logical :: test_OrbJast_psi_lapl  = .False.
 logical :: test_Jast    = .False.
 logical :: test_nucl_energy = .False.
 logical :: test_deriv_Jast  = .False.

 integer :: i,j,k,l,m
 double precision :: delta

 double precision,allocatable :: gradient(:,:), laplacian(:)
 double precision :: h

 real, allocatable :: P(:,:), P2(:,:)
 integer :: icount
!print *,  jast_U1_Simple
!print *,  jast_F_ee_value
!print *,  jast_F_eN_value
!print *,  jast_U1_density_lapl(1,2), jast_U1_density_lapl(1,1)
!print *,  jast_value
 PROVIDE ezfio_filename
 print *,  trim(ezfio_filename)
 
 allocate(P(ao_num,ao_num), P2(ao_num,ao_num))
 print *,  E_loc

 allocate(gradient(ao_num,4))

 double precision :: error
! do j=1,mo_num
!   print *,  j
BEGIN_TEMPLATE
  print *,  '$GRADIENT'

 gradient = 0.d0
 delta = 1.d-2
 print *,  ''
 do i=1,elec_num
  gradient(i,4) = 0.d0
  do l=1,3
   call update_qmckl_coord()
   TOUCH elec_coord
   h = elec_coord(i,l)
   gradient(i,l) = 0.d0
   gradient(i,4) +=              -205.d0/72.d0*( $VALUE $X)

   elec_coord(i,l) = h -4.*delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=                1.d0/280.d0*( $VALUE $X)
   gradient(i,4) +=                -1.d0/560.d0*( $VALUE $X)

   elec_coord(i,l) = h -3.*delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=               -4.d0/105.d0*( $VALUE $X)
   gradient(i,4) +=                 8.d0/315.d0*( $VALUE $X)

   elec_coord(i,l) = h -2.*delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=                1.d0/5.d0*( $VALUE $X)
   gradient(i,4) +=                -1.d0/5.d0*( $VALUE $X)

   elec_coord(i,l) = h -delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=               -4.d0/5.d0*( $VALUE $X)
   gradient(i,4) +=                 8.d0/5.d0*( $VALUE $X)

   elec_coord(i,l) = h +delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=                4.d0/5.d0*( $VALUE $X)
   gradient(i,4) +=                 8.d0/5.d0*( $VALUE $X)

   elec_coord(i,l) = h +2.*delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=               -1.d0/5.d0*( $VALUE $X)
   gradient(i,4) +=                -1.d0/5.d0*( $VALUE $X)

   elec_coord(i,l) = h +3.*delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=                4.d0/105.d0*( $VALUE $X)
   gradient(i,4) +=                 8.d0/315.d0*( $VALUE $X)

   elec_coord(i,l) = h +4.*delta
   call update_qmckl_coord()
   TOUCH elec_coord
   gradient(i,l) +=               -1.d0/280.d0*( $VALUE $X)
   gradient(i,4) +=                -1.d0/560.d0*( $VALUE $X)

   gradient(i,l) *= 1.d0/delta
   elec_coord(i,l) = h
   enddo
  gradient(i,4) *= 1.d0/(delta*delta)
 enddo
 call update_qmckl_coord()
 TOUCH elec_coord
 l = $L
 do i=1,elec_num
 !if ( ( $GRADIENT $Y/real(gradient(i,l)+1.e-12)>1.00005) .or. &
 !     ( $GRADIENT $Y/real(gradient(i,l)+1.e-12)<0.99995) ) then
   error = abs($GRADIENT $Y - real(gradient(i,l))+1.e-16)/($GRADIENT $Y+1.e-16)
   print '(I3,I3,I3,e16.6,e16.6,F16.3)', i,j,l, $GRADIENT $Y, real(gradient(i,l)), &
     error
 !endif
 enddo
 print *, ''

SUBST [ VALUE, X, GRADIENT, Y, L ]
 jast_value ; ; jast_grad_x ; (i) ; 1 ;;
 jast_value ; ; jast_grad_y ; (i) ; 2 ;;
 jast_value ; ; jast_grad_z ; (i) ; 3 ;;
 jast_value ; ; jast_lapl   ; (i) ; 4 ;;
END_TEMPLATE
! enddo



!mo_value_transp ; (j,i) ; mo_grad_transp_x ; (j,i) ; 1 ;;
!mo_value_transp ; (j,i) ; mo_grad_transp_y ; (j,i) ; 2 ;;
!mo_value_transp ; (j,i) ; mo_grad_transp_z ; (j,i) ; 3 ;;
!mo_value_transp ; (j,i) ; mo_lapl_transp   ; (j,i) ; 4 ;;

!ao_value_block ; (m) ; ao_grad_block_x ; (m) ; 1 ;;
!ao_value_block ; (m) ; ao_grad_block_y ; (m) ; 2 ;;
!ao_value_block ; (m) ; ao_grad_block_z ; (m) ; 3 ;;
!ao_value_block ; (m) ; ao_lapl_block   ; (m) ; 4 ;;


!  density_det         ; (i)  ; density_det_grad_x ; (i) ; 1 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_grad_x ; (i) ; 1 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_grad_y ; (i) ; 2 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_grad_z ; (i) ; 3 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_lapl   ; (i) ; 4 ;;
! jast_a ; (1,i) ; jast_a_lapl   ; (1,i) ; 4 ;;
! jast_fb_value ; (i) ; jast_fb_grad_x ; (i) ; 1 ;;
! jast_fb_value ; (i) ; jast_fb_grad_y ; (i) ; 2 ;;
! jast_fb_value ; (i) ; jast_fb_grad_z ; (i) ; 3 ;;
! jast_fb_value ; (i) ; jast_fb_lapl   ; (i) ; 4 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_grad_x ; (i) ; 1 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_grad_y ; (i) ; 2 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_grad_z ; (i) ; 3 ;;
! sum(jast_elec_Large_value) ; ; jast_elec_Large_lapl   ; (i) ; 4 ;;

!force_q(2,1) ; ; force_q_grad ; (i,1,2,1) ; 1 ;;
!force_q(2,1) ; ; force_q_grad ; (i,2,2,1) ; 2 ;;
!force_q(2,1) ; ; force_q_grad ; (i,3,2,1) ; 3 ;;
!force_q(2,1) ; ; force_q_lapl ; (i,2,1)   ; 4 ;;
!force_q(2,2) ; ; force_q_grad ; (i,1,2,2) ; 1 ;;
!force_q(2,2) ; ; force_q_grad ; (i,2,2,2) ; 2 ;;
!force_q(2,2) ; ; force_q_grad ; (i,3,2,2) ; 3 ;;
!force_q(2,2) ; ; force_q_lapl ; (i,2,2)   ; 4 ;;
!force_q(2,3) ; ; force_q_grad ; (i,1,2,3) ; 1 ;;
!force_q(2,3) ; ; force_q_grad ; (i,2,2,3) ; 2 ;;
!force_q(2,3) ; ; force_q_grad ; (i,3,2,3) ; 3 ;;
!force_q(2,3) ; ; force_q_lapl ; (i,2,3)   ; 4 ;;
!density_det         ; (i)  ; density_det_grad_x ; (i) ; 1 ;;
!density_det         ; (i)  ; density_det_grad_y ; (i) ; 2 ;;
!density_det         ; (i)  ; density_det_grad_z ; (i) ; 3 ;;
!density_det         ; (i)  ; density_det_lapl   ; (i) ; 4 ;;
!  sum(jast_F_ee_value);  ; jast_F_ee_grad_x ; (i) ; 1 ;;
!  sum(jast_F_ee_value);  ; jast_F_ee_grad_y ; (i) ; 2 ;;
!  sum(jast_F_ee_value);  ; jast_F_ee_grad_z ; (i) ; 3 ;;
!  sum(jast_F_ee_value);  ; jast_F_ee_lapl   ; (i) ; 4 ;;
! jast_F_eN_value;(i)  ; jast_F_eN_grad_x ; (i) ; 1 ;;
! jast_F_eN_value;(i)  ; jast_F_eN_grad_y ; (i) ; 2 ;;
! jast_F_eN_value;(i)  ; jast_F_eN_grad_z ; (i) ; 3 ;;
! jast_F_eN_value;(i)  ; jast_F_eN_lapl   ; (i) ; 4 ;;

! E_kin;  ; E_kin_grad ; (i,l) ;;
! E_loc;  ; E_loc_grad ; (i,l) ;;
! jast_grad(3,2);  ; jast_hess ; (j,l,3,2) ;;
! jast_lapl(3);  ; jast_lapl_grad ; (i,3,l) ;;
! E_loc;  ; E_loc_grad ; (i,l) ;;
! psidet_lapl_psidet_inv ; (j) ; psidet_lapl_psidet_inv_grad ; (i,j,l) ;;
! psi_lapl_psi_inv ; (j) ; psi_lapl_psi_inv_grad ; (i,j,l) ;;
! mo_lapl; (i,j) ; mo_lapl_grad; (i,j,l) ;;
! E_kin;  ; E_kin_grad ; (i,l) ;;
! E_pot;  ; E_pot_grad ; (i,l) ;;
! psidet_lapl ; (3) ; psidet_lapl_grad ; (i,3,l) ;;

! det_unique_grad ; (3,2) ; det_unique_hess ; (i,l,3,2,1) ;;
! psi_lapl ; (3) ; psi_lapl_grad ; (i,3,l) ;;
! E_kin;   ;
! E_kin_grad ; (i,l) ;;

! det_unique_beta_lapl_curr ; (8) ; det_unique_beta_lapl_grad_curr ; (i,8,l) ;;
! det_unique_beta_grad_curr ; (8,3) ; det_unique_beta_hess_curr ; (i,l,8,3) ;;
! ao_axis_lapl; (j,i) ;
! ao_axis_lapl_grad; (j,i,l) ;;



!print *,  ''
! do i=1,elec_num
!  print *,  jast_F_ee_lapl(i), laplacian(i)
! enddo

 return
end

!print *,  ''
!print *,  nucl_force_psi_grad(:)
!print *,  ''
!print *,  nucl_force_eloc_psi_grad(:)
!return

!allocate (force_q_grad_test (elec_num,nucl_num,3))
!do l=1,3
! do j=1,nucl_num
!  do i=1,elec_num
!    force_q_grad_test(i,j,l) = -6.d0*force_q(j,l)
!   enddo
!  enddo
!enddo
!do k=1,3
! do i=1,elec_num
!  elec_coord(i,k) += 0.05
!  TOUCH elec_coord
!  do l=1,3
!   do j=1,nucl_num
!    force_q_grad_test(i,j,l) += force_q(j,l)
!   enddo
!  enddo
!  elec_coord(i,k) -= 0.10
!  TOUCH elec_coord
!  do l=1,3
!   do j=1,nucl_num
!    force_q_grad_test(i,j,l) += force_q(j,l)
!   enddo
!  enddo
!  elec_coord(i,k) += 0.05
! enddo
!enddo
!do i=1,elec_num
!  do l=1,3
!   do j=1,nucl_num
!    force_q_grad_test(i,j,l) *= 1.d0/(0.05d0**2)
!   enddo
!  enddo
!enddo
!SOFT_TOUCH elec_coord

!do l=1,3 
! do k=1,nucl_num
!   do i=1,elec_num
!    print *, force_q_grad_test(i,k,l), force_q_lapl(i,k,l)
!   enddo
! enddo
!enddo

!stop 

!if (test_OrbJast_value) then
! print *,  'orb_jast_num', orb_jast_num
! print *,  'do_orb_jast', do_orb_jast
! print *,  'orb_jast_type', orb_jast_type
! print *,  'orb_jast_pen'
! do k=1,orb_jast_num
!  print *,  orb_jast_pen(:,k)
! enddo
! print *,  'orb_jast_map'
! do k=1,mo_num
!  print *,  orb_jast_map(k)
! enddo
! print *,  'orb_jast_a..b'
! do k=1,orb_jast_num
!  print *,  orb_jast_a(k), orb_jast_b(k)
! enddo
! print *,  ''
! double precision :: s
! s = 1.d0
! do k=1,elec_num
!  s = s * orb_jast_value(k,1)
! enddo
! print *, 'orb_jast_value', s
! print *,  'jast_a..b'
! print *,  jast_a_up_up,jast_a_up_dn,jast_b_up_up,jast_b_up_dn
! print *,  'jast_pen'
! print *,  jast_pen(:)
! print *, 'jast_value', jast_value
! print *,  ''
!endif

!if (test_OrbJast_grad) then
! print *, 'orb_jast_grad'
! double precision, allocatable :: orb_jast_grad_new(:,:,:)
! allocate(orb_jast_grad_new(elec_num,elec_num,3))
! do l=1,3
!  do j=1,elec_num
!   do i=1,elec_num
!    orb_jast_grad_new(i,j,l) = 0.d0
!   enddo
!  enddo
! enddo
! delta = 1.d-3
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) - delta
!   TOUCH elec_coord
!   do j=1,elec_num
!    orb_jast_grad_new(i,j,l) = orb_jast_grad_new(i,j,l) - orb_jast_value(j,1)
!   enddo
!   elec_coord(i,l) = elec_coord(i,l) + delta
!  enddo
! enddo
! TOUCH elec_coord
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) + delta
!   TOUCH elec_coord
!   do j=1,elec_num
!    orb_jast_grad_new(i,j,l) = orb_jast_grad_new(i,j,l) + orb_jast_value(j,1)
!   enddo
!   elec_coord(i,l) = elec_coord(i,l) - delta
!  enddo
! enddo
! TOUCH elec_coord
! 
! do l=1,3
!  do j=1,elec_num
!   do i=1,elec_num
!    orb_jast_grad_new(i,j,l)= orb_jast_grad_new(i,j,l)/(2.d0*delta)
!   enddo
!  enddo
! enddo

! print *,  '---'
! do j=1,elec_num
!  do i=1,elec_num
!   print *, orb_jast_grad_new(i,j,:)
!   print *, orb_jast_grad(i,j,:,1)
!   print *,  ''
!  enddo
! enddo
! print *,  '---'

! deallocate(orb_jast_grad_new)
!endif

!if (test_OrbJast_lapl) then
! print *, 'orb_jast_lapl'
! double precision, allocatable :: orb_jast_lapl_new(:,:)
! allocate(orb_jast_lapl_new(elec_num,elec_num))
! do j=1,elec_num
!  do i=1,elec_num
!   orb_jast_lapl_new(i,j) = -6.d0*orb_jast_value(j,1)
!  enddo
! enddo
! delta = 5.d-2
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) - delta
!   TOUCH elec_coord
!   do j=1,elec_num
!    orb_jast_lapl_new(i,j) = orb_jast_lapl_new(i,j) + orb_jast_value(j,1)
!   enddo
!   elec_coord(i,l) = elec_coord(i,l) + delta
!  enddo
! enddo
! TOUCH elec_coord
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) + delta
!   TOUCH elec_coord
!   do j=1,elec_num
!    orb_jast_lapl_new(i,j) = orb_jast_lapl_new(i,j) + orb_jast_value(j,1)
!   enddo
!   elec_coord(i,l) = elec_coord(i,l) - delta
!  enddo
! enddo
! TOUCH elec_coord
! 
! do j=1,elec_num
!  do i=1,elec_num
!   orb_jast_lapl_new(i,j)= orb_jast_lapl_new(i,j)/(delta**2)
!  enddo
! enddo

! print *,  '---'
! do j=1,elec_num
!  do i=1,elec_num
!   print *, orb_jast_lapl_new(i,j)
!   print *, orb_jast_lapl(i,j,1)
!   print *,  ''
!  enddo
! enddo
! print *,  '---'

! deallocate(orb_jast_lapl_new)
!endif


!if (test_OrbJast_psi_value) then
! print *,  ''
! do_orb_jast = .False.
! do_jast = .not.do_orb_jast
! TOUCH do_orb_jast do_jast
! print *, psi_value

! do_orb_jast = .True.
! do_jast = .not.do_orb_jast
! TOUCH do_orb_jast do_jast
! print *, psi_value
!endif

!if (test_OrbJast_psi_grad) then
! print *,  do_orb_jast, do_jast
! allocate(orb_jast_grad_new(elec_num,1,3))
! do l=1,3
!   do i=1,elec_num
!    orb_jast_grad_new(i,1,l) = 0.d0
!   enddo
! enddo
! delta = 1.d-4
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) - delta
!   TOUCH elec_coord
!   orb_jast_grad_new(i,1,l) = orb_jast_grad_new(i,1,l) - det_unique_alpha_value_curr
!   elec_coord(i,l) = elec_coord(i,l) + delta
!  enddo
! enddo
! TOUCH elec_coord
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) + delta
!   TOUCH elec_coord
!   orb_jast_grad_new(i,1,l) = orb_jast_grad_new(i,1,l) + det_unique_alpha_value_curr
!   elec_coord(i,l) = elec_coord(i,l) - delta
!  enddo
! enddo
! TOUCH elec_coord
! 
! do l=1,3
!  do i=1,elec_num
!   orb_jast_grad_new(i,1,l)= orb_jast_grad_new(i,1,l)/(2.d0*delta)
!  enddo
! enddo

! print *,  '---'
! do i=1,elec_num
!   print *, orb_jast_grad_new(i,1,:)
!   print *, det_unique_alpha_grad_curr(i,:)
!   print *,  ''
! enddo
! print *,  '---'

! deallocate(orb_jast_grad_new)
! print *,  ''
! do_orb_jast = .False.
! do_jast = .not.do_orb_jast
! TOUCH do_orb_jast do_jast
! do i=1,elec_num
!   print *, psi_grad(i,:)
! enddo
! print *,  E_loc

! do_orb_jast = .True.
! do_jast = .not.do_orb_jast
! TOUCH do_orb_jast do_jast
! print *, ''
! do i=1,elec_num
!   print *, psi_grad(i,:)
! enddo
! print *,  E_loc
! print *, ''
!endif

!if (test_OrbJast_psi_lapl) then
! allocate(orb_jast_lapl_new(elec_num,1))
! do i=1,elec_num
!    orb_jast_lapl_new(i,1) = -6.d0*det_unique_alpha_value_curr
! enddo
! delta = 5.d-2
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) - delta
!   TOUCH elec_coord
!   orb_jast_lapl_new(i,1) = orb_jast_lapl_new(i,1) + det_unique_alpha_value_curr
!   elec_coord(i,l) = elec_coord(i,l) + delta
!  enddo
! enddo
! TOUCH elec_coord
! do l=1,3
!  do i=1,elec_num
!   elec_coord(i,l) = elec_coord(i,l) + delta
!   TOUCH elec_coord
!   orb_jast_lapl_new(i,1) = orb_jast_lapl_new(i,1) + det_unique_alpha_value_curr
!   elec_coord(i,l) = elec_coord(i,l) - delta
!  enddo
! enddo
! TOUCH elec_coord
! 
! do i=1,elec_num
!   orb_jast_lapl_new(i,1)= orb_jast_lapl_new(i,1)/(delta**2)
! enddo

! print *,  '---'
! do i=1,elec_num
!   print *, orb_jast_lapl_new(i,1), det_unique_alpha_lapl_curr(i)
! enddo
! print *,  '---'

! deallocate(orb_jast_lapl_new)
! print *,  ''
! print *,  ''
! do_orb_jast = .False.
! do_jast = .not.do_orb_jast
! TOUCH do_orb_jast do_jast
! do i=1,elec_num
!   print *, psi_lapl(i)
! enddo

! do_orb_jast = .True.
! do_jast = .not.do_orb_jast
! TOUCH do_orb_jast do_jast
! print *, ''
! do i=1,elec_num
!   print *, psi_lapl(i)
! enddo
! print *, ''
!endif

!if (test_E_loc) then
! print *,  'E_loc'
! do k=1,walk_num

!  do l=1,3
!   do i=1,elec_num
!    elec_coord(i,l) = elec_coord_full(i,l,k)
!   enddo
!  enddo
!  TOUCH elec_coord

!  print *, E_loc
! enddo
!endif

!if (test_E_kin) then
! print *,  'E_kin'
! do k=1,walk_num

!  do l=1,3
!   do i=1,elec_num
!    elec_coord(i,l) = elec_coord_full(i,l,k)
!   enddo
!  enddo
!  TOUCH elec_coord

!  print *, E_kin
! enddo
!endif

!if (test_psi) then
! print *,  'psi'
! do k=1,walk_num

!  do l=1,3
!   do i=1,elec_num
!    elec_coord(i,l) = elec_coord_full(i,l,k)
!   enddo
!  enddo
!  TOUCH elec_coord

!  print *, psi_value
! enddo
!endif

!if (test_lapl) then
! print *,  'psi_lapl'
!     do i=1,elec_num
!       print *,  psi_lapl(i)
!     enddo
!     print *,  ''
!     double precision, allocatable :: psi_lapl_new(:)
!     allocate(psi_lapl_new(elec_num))
!     do i=1,elec_num
!      psi_lapl_new(i) = -6.d0*psi_value
!     enddo
!     delta = 5.d-2
!     do l=1,3
!      do i=1,elec_num
!       elec_coord(i,l) = elec_coord(i,l) - delta
!       TOUCH elec_coord
!       psi_lapl_new(i) = psi_lapl_new(i) + psi_value
!       elec_coord(i,l) = elec_coord(i,l) + delta
!      enddo
!     enddo
!     TOUCH elec_coord
!     do l=1,3
!      do i=1,elec_num
!       elec_coord(i,l) = elec_coord(i,l) + delta
!       TOUCH elec_coord
!       psi_lapl_new(i) = psi_lapl_new(i) + psi_value
!       elec_coord(i,l) = elec_coord(i,l) - delta
!      enddo
!     enddo
!     TOUCH elec_coord
!     
!     do i=1,elec_num
!       psi_lapl_new(i)= psi_lapl_new(i)/(delta**2)
!       print *,  psi_lapl_new(i)
!     enddo

!     deallocate(psi_lapl_new)
! endif

! if (test_Jast) then
!   print *,  'Test Jast'
!   print *,  '---------'

!   jast_type = t_Simple
!   TOUCH jast_type
!   print *,  'Simple : ', jast_value

!   jast_type = t_Spline
!   TOUCH jast_type
!   print *,  'Spline : ', jast_value

!   jast_type = t_Simple
!   TOUCH jast_type
!   print *,  'Simple : '
!   do l=1,3
!    do i=1,elec_num
!     print *,  i,l, jast_grad(i,l)
!    enddo
!   enddo
!   jast_type = t_Spline
!   TOUCH jast_type
!   print *,  'Spline : '
!   do l=1,3
!    do i=1,elec_num
!     print *,  i,l, jast_grad(i,l)
!    enddo
!   enddo

!   jast_type = t_Simple
!   TOUCH jast_type
!   print *,  'Simple : '
!   do i=1,elec_num
!    print *,  i, jast_lapl(i)
!   enddo
!   jast_type = t_Spline
!   TOUCH jast_type
!   print *,  'Spline : '
!   do i=1,elec_num
!     print *,  i, jast_lapl(i)
!   enddo

!   jast_type = t_Simple
!   TOUCH jast_type
!   print *,  'E_simple : ', E_loc

!   jast_type = t_Spline
!   TOUCH jast_type
!   print *,  'E_spline : ', E_loc

!   print *,  ''
! endif

! if (test_nucl_energy) then
!   print *, E_nucl
! endif

! print *,  '----'

! if (test_nucl_energy) then
!   print *, E_nucl
! endif

! print *,  '----'
! if (test_deriv_Jast) then
!   print *,  d_ene_jast_een_e_b
!   double precision :: deriv
!   delta = 1.d-3
!   deriv = 0.d0
!   jast_een_e_b -= delta
!   TOUCH jast_een_e_b
!   deriv -= psi_value
!   jast_een_e_b += 2.d0*delta
!   TOUCH jast_een_e_b
!   deriv += psi_value
!   deriv = deriv / (2.d0*delta)
!   jast_een_e_b -= delta
!   TOUCH jast_een_e_b
!   print *,  deriv/psi_value
! endif
! print *,  '----'

