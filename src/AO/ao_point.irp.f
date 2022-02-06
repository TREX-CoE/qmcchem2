BEGIN_PROVIDER [ real, ao_axis_power_p, (-2:ao_power_max,3,nucl_num) ]
  implicit none
  BEGIN_DOC
! Evaluation of power of x, y, z at the current point for each
! nucleus. Negative power -> 0.
  END_DOC
  
  integer                        :: i,k,l
  do i=1,nucl_num
    do l=1,3
      ao_axis_power_p(-2,l,i) = 0.
      ao_axis_power_p(-1,l,i) = 0.
      ao_axis_power_p(0,l,i) = 0.
      ao_axis_power_p(0,l,i) = 1.
      do k=1,ao_power_max_nucl(i,l)
        ao_axis_power_p(k,l,i) = point_nucl_dist_vec(i,l)*ao_axis_power_p(k-1,l,i)
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_axis_p, (ao_num) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Cartesian polynomial part of the atomic orbitals.
  END_DOC
  integer                        :: i
  
  do i=1,ao_num
    ao_axis_p(i)                                                     &
        = ao_axis_power_p( ao_power_transp(1,i) , 1 , ao_nucl(i) )   &
        * ao_axis_power_p( ao_power_transp(2,i) , 2 , ao_nucl(i) )   &
        * ao_axis_power_p( ao_power_transp(3,i) , 3 , ao_nucl(i) )
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_axis_grad_p, (ao_num,3) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Gradients of the cartesian polynomial part of the atomic orbitals.
  END_DOC
  integer                        :: i, l
  real                           :: real_of_int(-1:10)
  data real_of_int /0.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
  
  do i=1,ao_num
    ao_axis_grad_p(i,1) = real_of_int(ao_power_transp(1,i))          &
        * ao_axis_power_p( ao_power_transp(1,i)-1, 1 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(2,i)  , 2 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(3,i)  , 3 , ao_nucl(i) )
  enddo
  
  do i=1,ao_num
    ao_axis_grad_p(i,2) = real_of_int(ao_power_transp(2,i))          &
        * ao_axis_power_p( ao_power_transp(1,i)  , 1 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(2,i)-1, 2 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(3,i)  , 3 , ao_nucl(i) )
  enddo
  
  do i=1,ao_num
    ao_axis_grad_p(i,3) = real_of_int(ao_power_transp(3,i))          &
        * ao_axis_power_p( ao_power_transp(1,i)  , 1 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(2,i)  , 2 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(3,i)-1, 3 , ao_nucl(i) )
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_axis_lapl_p, (ao_num) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Laplacian of the cartesian atomic orbitals 
  END_DOC
  integer                        :: i, j, l
  
  do i=1,ao_num
    real                           :: real_of_int(-2:10)
    data real_of_int /0.,0.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
    
    ao_axis_lapl_p(i)                                                &
        = real_of_int(ao_power_transp(1,i))                          &
        * real_of_int(ao_power_transp(1,i)-1)                        &
        * ao_axis_power_p( ao_power_transp(1,i)-2, 1 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(2,i)  , 2 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(3,i)  , 3 , ao_nucl(i) )  &
        + real_of_int(ao_power_transp(2,i))                          &
        * real_of_int(ao_power_transp(2,i)-1)                        &
        * ao_axis_power_p( ao_power_transp(1,i)  , 1 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(2,i)-2, 2 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(3,i)  , 3 , ao_nucl(i) )  &
        + real_of_int(ao_power_transp(3,i))                          &
        * real_of_int(ao_power_transp(3,i)-1)                        &
        * ao_axis_power_p( ao_power_transp(1,i)  , 1 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(2,i)  , 2 , ao_nucl(i) )  &
        * ao_axis_power_p( ao_power_transp(3,i)-2, 3 , ao_nucl(i) )
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_value_p, (ao_num) ]
  implicit none
  BEGIN_DOC
! Values of the atomic orbitals
  END_DOC
  integer :: i
  do i=1,ao_num
    ao_value_p(i) = ao_oneD_p(i) * ao_axis_p(i)
  enddo
END_PROVIDER


BEGIN_PROVIDER [ real, ao_grad_p, (ao_num,3) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Gradients of the atomic orbitals 
  END_DOC

  integer                        :: i,l
  do l=1,3
    do i=1,ao_num
      ao_grad_p(i,l) = ao_oneD_p(i) * ao_axis_grad_p(i,l)
    enddo
    do i=1,ao_num
      ao_grad_p(i,l) = ao_grad_p(i,l) + ao_oneD_grad_p(i,l) * ao_axis_p(i)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_lapl_p, (ao_num) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Laplacian of the atomic orbitals 
  END_DOC
  integer                        :: i,l
  do i=1,ao_num
    ao_lapl_p(i) = ao_oneD_p(i) * ao_axis_lapl_p(i)
  enddo
  do i=1,ao_num
    ao_lapl_p(i) = ao_lapl_p(i) + ao_oneD_lapl_p(i) * ao_axis_p(i)
  enddo
  do l=1,3
    do i=1,ao_num
      ao_lapl_p(i) = ao_lapl_p(i) + 2.*ao_oneD_grad_p(i,l) * ao_axis_grad_p(i,l)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_oneD_prim_p, (ao_num,ao_prim_num_max) ]
  implicit none
  include '../types.F'

  BEGIN_DOC
! Exponentials of the primitive AOs
  END_DOC
  integer                        :: i,  k
  real                           :: r2, rtemp
  
  ! Compute alpha*r or alpha*r^2
  do i=1,ao_num
    r2 = point_nucl_dist(ao_nucl(i))*point_nucl_dist(ao_nucl(i))
    do k=1,ao_prim_num_max
      ao_oneD_prim_p(i,k) = r2
    enddo
  enddo
  
  ! Compute exp(-alpha*r) or exp(-alpha*r^2)
  do i=1,ao_num
    do k=1,ao_prim_num(i)
      ao_oneD_prim_p(i,k) = exp(-ao_oneD_prim_p(i,k)*ao_expo(i,k))
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_oneD_p, (ao_num) ]
  implicit none
  include '../types.F'

  BEGIN_DOC
! One-dimensional component of the AOs
  END_DOC
  
  integer                        :: i, k
  
  do i=1,ao_num
    ao_oneD_p(i) = 0.
  enddo
  do k=1,ao_prim_num_max
    do i=1,ao_num
      ao_oneD_p(i) = ao_oneD_p(i) + ao_coef(i,k)*ao_oneD_prim_p(i,k)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_oneD_prim_grad_p, (ao_num,ao_prim_num_max,3) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Gradients of the one-dimensional component of the primitive AOs
  END_DOC
  integer                        :: i, k, l
  do l=1,3
    do k=1,ao_prim_num_max
      do i=1,ao_num
        ao_oneD_prim_grad_p(i,k,l) = -2.*point_nucl_dist_vec(ao_nucl(i),l)*ao_expo(i,k)*ao_oneD_prim_p(i,k)
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_oneD_grad_p, (ao_num,3) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Gradients of the one-dimensional component of the AOs
  END_DOC
  integer                        :: i, k, l
  do l=1,3
    do i=1,ao_num
      ao_oneD_grad_p(i,l) = 0.
    enddo
    do k=1,ao_prim_num_max
      do i=1,ao_num
        ao_oneD_grad_p(i,l) = ao_oneD_grad_p(i,l) + ao_coef(i,k)*ao_oneD_prim_grad_p(i,k,l)
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_oneD_prim_lapl_p, (ao_num,ao_prim_num_max) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Laplacian of the one-dimensional component of the primitive AOs
  END_DOC
  integer                        :: i, k
  do k=1,ao_prim_num_max
    do i=1,ao_num
      ao_oneD_prim_lapl_p(i,k) = ao_oneD_prim_p(i,k) * ao_expo(i,k) *&
          ( 4.*ao_expo(i,k)*point_nucl_dist(ao_nucl(i))**2 - 6. )
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ real, ao_oneD_lapl_p, (ao_num) ]
  implicit none
  include '../types.F'
  BEGIN_DOC
! Laplacian of the one-dimensional component of the AOs
  END_DOC
  
  integer                        :: i, k
  
  do i=1,ao_num
    ao_oneD_lapl_p(i) = 0.
  enddo
  do k=1,ao_prim_num_max
    do i=1,ao_num
      ao_oneD_lapl_p(i) = ao_oneD_lapl_p(i) + ao_coef(i,k)*ao_oneD_prim_lapl_p(i,k)
    enddo
  enddo

END_PROVIDER


