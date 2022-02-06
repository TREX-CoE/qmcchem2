subroutine pow_l(r,a,x1,x2,x3)
  implicit none
  BEGIN_DOC
! Fast calculation of powers for AO calculations
  END_DOC
  real, intent(in) :: r
  integer, intent(in) :: a
  real, intent(out) :: x1,x2,x3
  if (a==0) then
    x1 = 1.
    x2 = 0.
    x3 = 0.
     return
  else if (a==1) then
    x1 = r
    x2 = 1.
    x3 = 0.
     return
  else if (a==2) then
    x1 = r*r
    x2 = r
    x3 = 1.
    return
  endif

  select case (a)
   case (3)
     x2 = r*r
     x3 = r
     x1 = x2*r
     return
   case (4)
     x2 = r*r*r
     x3 = r*r
     x1 = x3*x3
     return
   case (5:)
     x3 = r**(a-2)
     x2 = x3*r
     x1 = x2*r
     return
   case (:-1)
     x1 = 0.
     x2 = 0.
     x3 = 0.
     return
  end select
end


 BEGIN_PROVIDER [ real, ao_axis_block, (ao_block_num_8) ]
&BEGIN_PROVIDER [ real, ao_axis_grad_block_x, (ao_block_num_8) ]
&BEGIN_PROVIDER [ real, ao_axis_grad_block_y, (ao_block_num_8) ]
&BEGIN_PROVIDER [ real, ao_axis_grad_block_z, (ao_block_num_8) ]
&BEGIN_PROVIDER [ real, ao_axis_lapl_block, (ao_block_num_8) ]
  implicit none
  include '../types.F'

  BEGIN_DOC
! Cartesian polynomial part of the atomic orbitals. (blocked)

! Gradients of the cartesian polynomial part of the atomic orbitals. (blocked)

! Laplacian of the cartesian atomic orbitals (blocked)
  END_DOC


  integer :: i, j, k, l, idx
  real, save :: real_of_int(-1:15)
  data real_of_int /0.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15./
  !DIR$ ATTRIBUTES ALIGN : 64 :: real_of_int

  if (ao_block_num_8 == 0) then
    return
  endif

  !DIR$ VECTOR ALIGNED
  do idx=1,ao_oneD_prim_non_zero_idx(0)
   i=ao_oneD_prim_non_zero_idx(idx)

   real :: p10, p11, p12   ! p10 contains x**(ao_power(i,1))
   real :: p20, p21, p22   ! p11 contains x**(ao_power(i,1)-1)
   real :: p30, p31, p32   ! p12 contains x**(ao_power(i,1)-2)
   real :: p012, p013, p023

   integer :: inucl
   inucl = ao_nucl(i)

   real :: x, y, z
   x = nucl_elec_dist_vec(1,inucl,ao_elec)
   y = nucl_elec_dist_vec(2,inucl,ao_elec)
   z = nucl_elec_dist_vec(3,inucl,ao_elec)

   integer :: pow1, pow2, pow3

   pow1 = ao_power_transp(1,i)
   pow2 = ao_power_transp(2,i)
   pow3 = ao_power_transp(3,i)

   !DIR$ FORCEINLINE
   call pow_l(x,pow1,p10, p11, p12)
   !DIR$ FORCEINLINE
   call pow_l(y,pow2,p20, p21, p22)
   !DIR$ FORCEINLINE
   call pow_l(z,pow3,p30, p31, p32)

   p012 = real_of_int(pow3) * p10*p20
   p023 = p20*p30
   p013 = real_of_int(pow2) * p10*p30

   ao_axis_block(idx) = p023 * p10
   p023 = real_of_int(pow1) * p023

   ao_axis_grad_block_x(idx) =  p023 * p11 
   ao_axis_grad_block_y(idx) =  p013 * p21
   ao_axis_grad_block_z(idx) =  p012 * p31
   ao_axis_lapl_block(idx) = real_of_int(pow1-1) * p023 * p12 & 
                           + real_of_int(pow2-1) * p013 * p22 &
                           + real_of_int(pow3-1) * p012 * p32



  enddo


END_PROVIDER



