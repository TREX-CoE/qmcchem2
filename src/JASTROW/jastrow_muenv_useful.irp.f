! Mu Jastrow x envelop
! ---------------------

! ---

BEGIN_PROVIDER [ integer, List_all_comb_b2_size]

  implicit none

  List_all_comb_b2_size = 2**nucl_num

END_PROVIDER

! ---

BEGIN_PROVIDER [ integer, List_all_comb_b2, (nucl_num, List_all_comb_b2_size)]

  implicit none
  integer :: i, j

  if(nucl_num .gt. 32) then
    print *, ' nucl_num = ', nucl_num, '> 32'
    stop
  endif

  List_all_comb_b2 = 0

  do i = 0, List_all_comb_b2_size-1
    do j = 0, nucl_num-1
      if (btest(i,j)) then
        List_all_comb_b2(j+1,i+1) = 1
      endif
    enddo
  enddo

END_PROVIDER

! ---

subroutine j_elec_Muenv(r1, r2, je)

  BEGIN_DOC  
  !
  ! J(i,j) = 0.5 [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ] 
  !        x v(riA) x v(rjA)
  !
  END_DOC

  include '../constants.F'

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: je
  double precision              :: mu_rij
  double precision              :: dx, dy, dz, rij, u_ij, vi, vj

  dx     = r1(1) - r2(1)
  dy     = r1(2) - r2(2)
  dz     = r1(3) - r2(3)
  rij    = dsqrt(dx*dx + dy*dy + dz*dz)
  mu_rij = mu_erf * rij
  u_ij   = 0.5d0 * (rij * (1.d0 - derf(mu_rij)) - dexp(-mu_rij*mu_rij)/(dsqpi*mu_erf))

  call j_elec_env(r1, vi)
  call j_elec_env(r2, vj)

  je = u_ij * vi * vj

  return
end subroutine j_elec_Muenv

! ---

subroutine j_elec_env(r, j1b)

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: j1b
  integer                       :: iA
  double precision              :: a, riA, r2
  double precision              :: dx, dy, dz

  PROVIDE j1b_type

  if((j1b_type .eq. 2) .or. (j1b_type .eq. 102)) then

    j1b = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = j1b_pen(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      riA = dsqrt(dx*dx + dy*dy + dz*dz)
      j1b = j1b - dexp(-a*riA)
    enddo

  elseif((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    j1b = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = j1b_pen(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      riA = dsqrt(dx*dx + dy*dy + dz*dz)
      j1b = j1b * (1.d0 - dexp(-a*riA*riA))
    enddo

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    j1b = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = j1b_pen(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      r2  = dx*dx + dy*dy + dz*dz
      j1b = j1b - dexp(-a*r2)
    enddo

  elseif((j1b_type .eq. 5) .or. (j1b_type .eq. 105)) then

    j1b = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = j1b_pen(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      r2  = dx*dx + dy*dy + dz*dz
      j1b = j1b - dexp(-a*r2*r2)
    enddo

  else

    print*, 'j1_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end subroutine j_elec_env

! ---

