BEGIN_PROVIDER [ double precision, pseudo_v_k , (nucl_num,pseudo_klocmax) ]
  implicit none
  BEGIN_DOC
! V_k
  END_DOC
  call get_pseudo_pseudo_v_k(pseudo_v_k)

END_PROVIDER

BEGIN_PROVIDER [ double precision, pseudo_v_kl , (nucl_num,pseudo_kmax,0:pseudo_lmax) ]
  implicit none
  BEGIN_DOC
! V_Kl
  END_DOC
  call get_pseudo_pseudo_v_kl(pseudo_v_kl)

END_PROVIDER

BEGIN_PROVIDER [ integer, pseudo_kmax  ]
  implicit none
  BEGIN_DOC
! kmax
  END_DOC
  call get_pseudo_pseudo_kmax(pseudo_kmax)

END_PROVIDER

BEGIN_PROVIDER [ double precision, pseudo_dz_kl , (nucl_num,pseudo_kmax,0:pseudo_lmax) ]
  implicit none
  BEGIN_DOC
! Exponents in the non-local part of the pseudo-potential
  END_DOC
  call get_pseudo_pseudo_dz_kl(pseudo_dz_kl)

END_PROVIDER

BEGIN_PROVIDER [ integer, pseudo_klocmax  ]
  implicit none
  BEGIN_DOC
! klocmax
  END_DOC
  call get_pseudo_pseudo_klocmax(pseudo_klocmax)

END_PROVIDER

BEGIN_PROVIDER [ integer, pseudo_lmax  ]
  implicit none
  BEGIN_DOC
! Max value of l in the pseudo-potential
  END_DOC
  call get_pseudo_pseudo_lmax(pseudo_lmax)

END_PROVIDER

BEGIN_PROVIDER [ integer, pseudo_grid_size  ]
  implicit none
  BEGIN_DOC
! Size of the QMC grid (number of points)
  END_DOC
  call get_pseudo_pseudo_grid_size(pseudo_grid_size)

END_PROVIDER

BEGIN_PROVIDER [ double precision, pseudo_grid_rmax  ]
  implicit none
  BEGIN_DOC
! Size of the QMC grid (max distance)
  END_DOC
  call get_pseudo_pseudo_grid_rmax(pseudo_grid_rmax)

END_PROVIDER

 BEGIN_PROVIDER [ integer, pseudo_non_loc_dim ]
&BEGIN_PROVIDER [ integer, pseudo_non_loc_dim_8 ]
&BEGIN_PROVIDER [ integer, pseudo_non_loc_dim_count, (nucl_num) ]
 implicit none
 BEGIN_DOC
 ! Dimension of of pseudo_non_local arrays
 END_DOC
 pseudo_non_loc_dim = 0
 pseudo_non_loc_dim_count = 0
 integer :: k,l,m
 do k=1,nucl_num
   do l=0,pseudo_lmax
     pseudo_non_loc_dim_count(k) += 2*l+1
   enddo
 enddo
 pseudo_non_loc_dim = sum(pseudo_non_loc_dim_count)
 integer, external :: mod_align
 pseudo_non_loc_dim_8 = mod_align(pseudo_non_loc_dim)
END_PROVIDER

BEGIN_PROVIDER [ integer, pseudo_n_kl , (nucl_num,pseudo_kmax,0:pseudo_lmax) ]
  implicit none
  BEGIN_DOC
! n_kl
  END_DOC
  call get_pseudo_pseudo_n_kl(pseudo_n_kl)

END_PROVIDER

BEGIN_PROVIDER [ double precision, pseudo_dz_k , (nucl_num,pseudo_klocmax) ]
  implicit none
  BEGIN_DOC
! dz_k
  END_DOC
  call get_pseudo_pseudo_dz_k(pseudo_dz_k)

END_PROVIDER

BEGIN_PROVIDER [ integer, pseudo_n_k , (nucl_num,pseudo_klocmax) ]
  implicit none
  BEGIN_DOC
! n_k
  END_DOC
  call get_pseudo_pseudo_n_k(pseudo_n_k)

END_PROVIDER

BEGIN_PROVIDER [ logical, do_pseudo  ]
  implicit none
  BEGIN_DOC
! Using pseudo potential integral of not
  END_DOC
  call get_pseudo_do_pseudo(do_pseudo)

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_pseudo_grid, (ao_num, -pseudo_lmax:pseudo_lmax, 0:pseudo_lmax, nucl_num, pseudo_grid_size) ]
 implicit none
 BEGIN_DOC
 ! Pseudopotential grid points
 END_DOC
 call get_pseudo_ao_pseudo_grid(ao_pseudo_grid)
END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_pseudo_grid, (ao_num, -pseudo_lmax:pseudo_lmax, 0:pseudo_lmax, nucl_num, pseudo_grid_size) ]
 implicit none
 BEGIN_DOC
 ! Pseudopotential grid points
 END_DOC
 call get_pseudo_mo_pseudo_grid(mo_pseudo_grid)
END_PROVIDER


BEGIN_PROVIDER [ double precision, mo_pseudo_grid_scaled, (pseudo_non_loc_dim_8,ao_num,pseudo_grid_size) ]
 implicit none
 BEGIN_DOC
 ! Pseudopotential grid points
 END_DOC
 integer                        :: i,k,l,m,kk,n
 double precision               :: c
 c = 1.d0/mo_scale
 do n=1,pseudo_grid_size
   do i=1,ao_num
     kk = 0
     do k=1,nucl_num
       do l=0,pseudo_lmax
         do m=-l,l
           kk = kk+1
           mo_pseudo_grid_scaled(kk,i,n) = c * mo_pseudo_grid(i,m,l,k,n)
         enddo
       enddo
     enddo
   enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, v_pseudo_local, (elec_num) ]
 implicit none
 BEGIN_DOC
 ! Local component of the pseudo-potential
 END_DOC
 integer :: i,k,l
 double precision :: r, rn, alpha
 do l=1,elec_num
   v_pseudo_local(l) = 0.d0
   do k=1,pseudo_klocmax
     do i=1,nucl_num
       r = nucl_elec_dist(i,l)
       alpha = pseudo_dz_k(i,k)*r*r
       if (alpha > 20.d0) then
         cycle
       endif
       select case (pseudo_n_k(i,k))
         case (-3)
           rn = nucl_elec_dist_inv(i,l)*nucl_elec_dist_inv(i,l)*nucl_elec_dist_inv(i,l)
         case (-2)
           rn = nucl_elec_dist_inv(i,l)*nucl_elec_dist_inv(i,l)
         case (-1)
           rn = nucl_elec_dist_inv(i,l)
         case (0)
           rn = 1.d0
         case (1)
           rn = r
         case (2)
           rn = r*r
         case (3)
           rn = r*r*r
         case default
           rn = r**pseudo_n_k(i,k)
       end select
       v_pseudo_local(l) += pseudo_v_k(i,k) * exp(-alpha) * rn
     enddo
   enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, v_pseudo_non_local, (pseudo_non_loc_dim_8,elec_num) ]
 implicit none
 BEGIN_DOC
 ! Non-Local component of the pseudo-potential
 END_DOC
 integer                        :: i,k,j,l
 integer                        :: kk,m
 double precision               :: r, rn, r2, alpha
 double precision               :: tmp(0:pseudo_lmax)
 v_pseudo_non_local = 0.d0
 do j=1,elec_num

   kk = 0
   do i=1,nucl_num
     r = nucl_elec_dist(i,j)
     r2 = r*r
     tmp(0:pseudo_lmax) = 0.d0
     do k=1,pseudo_kmax
       do l=0,pseudo_lmax
         alpha = pseudo_dz_kl(i,k,l)*r2
         if (alpha > 20.d0) then
           cycle
         endif
         select case (pseudo_n_kl(i,k,l))
           case (0)
             rn = 1.d0
           case (1)
             rn = r
           case (-1)
             rn = nucl_elec_dist_inv(i,j)
           case (2)
             rn = r2
           case (-2)
             rn = nucl_elec_dist_inv(i,j)*nucl_elec_dist_inv(i,j)
           case (3)
             rn = r2*r
           case (-3)
             rn = nucl_elec_dist_inv(i,j)*nucl_elec_dist_inv(i,j)*nucl_elec_dist_inv(i,j)
           case default
             rn = r**pseudo_n_kl(i,k,l)
         end select
         tmp(l) += pseudo_v_kl(i,k,l) * dexp(-alpha) * rn
       enddo
     enddo
     do l=0,pseudo_lmax
       do m=-l,l
         kk += 1
         v_pseudo_non_local(kk,j) = tmp(l)
       enddo
     enddo
   enddo

 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, pseudo_mo_term, (mo_num,elec_num) ]
 implicit none
 BEGIN_DOC
! \sum < Ylm | MO > x Ylm(r) x V_nl(r)
 END_DOC
 integer :: ii,i,j,l,m,k,n,kk
 double precision :: r, w0, w1, w2, ndr
 double precision :: tmp(pseudo_non_loc_dim_8,mo_num)
 double precision, save :: dr_inv, dr
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: tmp
 integer, save                  :: ifirst = 0

 if (ifirst == 0) then
    pseudo_mo_term = 0.d0
    ifirst = 1
    dr_inv = dble(pseudo_grid_size)/pseudo_grid_rmax
    dr     = pseudo_grid_rmax/dble(pseudo_grid_size)
 endif

 PROVIDE pseudo_ylm present_mos mo_pseudo_grid_scaled pseudo_non_loc_dim_count
 do j=1,elec_num
   tmp = 0.d0
   kk=0
   do k=1,nucl_num
     r = nucl_elec_dist(k,j)
     n = 1 + int(0.5+r*dr_inv)
     if (n<=pseudo_grid_size) then
       if (n == pseudo_grid_size) then
         n = n-1
       else if (n == 1) then
         n = 2
       endif

       ndr = dble(n-1)*dr
       w0 = -(r-ndr+dr)*dr_inv
       w1 = -(r-ndr)*dr_inv*0.5d0
       w2 = -(r-ndr-dr)*dr_inv

       do ii=1,num_present_mos
         i = present_mos(ii)
         !DIR$ LOOP COUNT(4)
         do l=kk+1,kk+pseudo_non_loc_dim_count(k)
           tmp(l,i) = ( ( mo_pseudo_grid_scaled (l,i,n-1) * w1 -        &
                          mo_pseudo_grid_scaled (l,i,n  ) * w0 ) * w2 + &
                        mo_pseudo_grid_scaled   (l,i,n+1) * w1 * w0 )   &
                      * pseudo_ylm(l,j)
         enddo
       enddo
     endif
     kk = kk+pseudo_non_loc_dim_count(k)
   enddo
   do ii=1,num_present_mos
     i = present_mos(ii)
     pseudo_mo_term(i,j) = 0.d0
     do kk=1,pseudo_non_loc_dim
       pseudo_mo_term(i,j) = pseudo_mo_term(i,j) + tmp(kk,i) * v_pseudo_non_local(kk,j)
     enddo
   enddo
 enddo
END_PROVIDER


double precision function ylm(l,m,x_,y_,z_,r_inv)
  implicit none
  include 'constants.F'
  BEGIN_DOC
! Y(l,m) evaluated at x,y,z
  END_DOC
  integer, intent(in)            :: l,m
  real, intent(in)               :: x_,y_,z_,r_inv
  double precision               :: x, y, z
  x = x_
  y = y_
  z = z_

  ylm = 0.5d0 * one_over_sqpi
  select case(l)

    case(0)
      continue

    case(1)
      ylm = ylm * sq3 * r_inv
      select case (m)
        case (-1)
          ylm = ylm * y
        case (0)
          ylm = ylm * z
        case (1)
          ylm = ylm * x
      end select

    case(2)
      ylm = ylm * r_inv * r_inv * dsqrt(5.d0)
      select case (m)
        case(-2)
          ylm =  ylm * sq3 * x * y
        case(-1)
          ylm =  ylm * sq3 * y * z
        case(0)
!         ylm =  ylm * 0.5d0 * (2.d0*z*z - x*x - y*y)
!          ylm =  ylm * 0.5d0 * (z*z - x*x + z*z - y*y)
          ylm =  ylm * 0.5d0 * ( (z+x)*(z-x) + (z+y)*(z-y) )
        case(1)
          ylm =  ylm * sq3 * z * x
        case(2)
!          ylm =  ylm * 0.5d0 * sq3 * (x*x - y*y)
          ylm =  ylm * 0.5d0 * sq3 * (x+y)*(x-y)
      end select

    case(3)
      ylm = ylm * r_inv * r_inv * r_inv * dsqrt(7.d0)
      select case (m)
        case(-3)
!          ylm =  ylm * 0.25d0 * dsqrt(10.d0) * (3.d0*x*x - y*y)
          ylm =  ylm * 0.25d0 * dsqrt(10.d0) * (2.d0*x*x + (x+y)*(x-y))
        case(-2)
          ylm =  ylm * dsqrt(15.d0) * x*y*z
        case(-1)
!          ylm =  ylm * 0.25d0*dsqrt(6.d0) * y * (4.d0*z*z - x*x -y*y)
          ylm =  ylm * 0.25d0*dsqrt(6.d0) * y * (2.d0*z*z + (z+x)*(z-x) + (z+y)*(z-y))
        case(0)
!          ylm =  ylm * 0.5d0 * z * (2.d0*z*z - 3.d0*(x*x + y*y))
!          ylm =  ylm * 0.5d0 * z * ( (z+x)*(z-x) + (z+y)*(z-y) - 2.d0*(x*x + y*y))
          ylm =  ylm * 0.5d0 * z * ( (z+x)*(z-x) + (z+y)*(z-y) - 2.d0*(x+y)*(x-y) - 4.d0*y*y)
        case(1)
!          ylm =  ylm * 0.25d0*dsqrt(6.d0) * x * (4.d0*z*z - x*x -y*y)
          ylm =  ylm * 0.25d0*dsqrt(6.d0) * x * (2.d0*z*z + (z-x)*(z+x) + (z-y)*(z+y))
        case(2)
!          ylm =  ylm * 0.5d0*dsqrt(15.d0) * z * (x*x -y*y)
          ylm =  ylm * 0.5d0*dsqrt(15.d0) * z * (x+y)*(x-y)
        case(3)
!          ylm =  ylm * 0.25d0*dsqrt(10.d0) * x * (x*x - 3.d0 * y*y)
          ylm =  ylm * 0.25d0*dsqrt(10.d0) * x * ( (x+y)*(x-y) - 2.d0*y*y)
      end select


    case default
      print *,  'l=', l
      stop 'problem in Ylm of pseudo : Ylm not implemented (pseudo.irp.f)'

  end select

end



BEGIN_PROVIDER [ double precision, pseudo_ylm, (pseudo_non_loc_dim_8,elec_num) ]
  implicit none
  BEGIN_DOC
  ! Y(l,m) evaluated for every electron position centered on every nuclei
  END_DOC

  integer                        :: i,j,l,m,kk
  double precision, external     :: ylm
  integer, save                  :: ifirst = 0

  if (ifirst == 0) then
    pseudo_ylm = 0.d0
    ifirst = 1
  endif

  do j=1,elec_num
    kk = 0
    do i=1,nucl_num
      do l=0,pseudo_lmax
        do m=-l,l
          kk = kk+1
          pseudo_ylm(kk,j) = ylm(l,m,                                &
              nucl_elec_dist_vec(1,i,j),                             &
              nucl_elec_dist_vec(2,i,j),                             &
              nucl_elec_dist_vec(3,i,j),                             &
              nucl_elec_dist_inv(i,j))
        enddo
      enddo
    enddo
  enddo
END_PROVIDER

