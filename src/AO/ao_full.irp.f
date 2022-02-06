 BEGIN_PROVIDER [ real, ao_value_block, (ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_grad_block_x, (ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_grad_block_y, (ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_grad_block_z, (ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_lapl_block, (ao_num_8) ]
&BEGIN_PROVIDER [ integer, ao_value_non_zero_idx, ((-simd_sp+1):ao_num_8) ]

  implicit none
  BEGIN_DOC
! Values of the atomic orbitals (blocked)

! Gradients of the atomic orbitals  (blocked)

! Laplacian of the atomic orbitals  (blocked)
  END_DOC

  if (ao_block_num == 0) then
    return
  endif

  integer :: j,idx,k, kmax
  ao_value_non_zero_idx(0) = ao_oneD_prim_non_zero_idx(0)
  kmax = ao_oneD_prim_non_zero_idx(0)
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do idx=1,kmax
   ao_value_non_zero_idx(idx) = ao_oneD_prim_non_zero_idx(idx)
   ao_value_block(idx) = ao_oneD_block(idx) * ao_axis_block(idx)
   ao_grad_block_x(idx) = ao_oneD_block(idx) * ao_axis_grad_block_x(idx) + ao_oneD_grad_block_x(idx) * ao_axis_block(idx) 
   ao_grad_block_y(idx) = ao_oneD_block(idx) * ao_axis_grad_block_y(idx) + ao_oneD_grad_block_y(idx) * ao_axis_block(idx) 
   ao_grad_block_z(idx) = ao_oneD_block(idx) * ao_axis_grad_block_z(idx) + ao_oneD_grad_block_z(idx) * ao_axis_block(idx) 
   ao_lapl_block(idx) = ao_oneD_block(idx) * ao_axis_lapl_block(idx) + ao_oneD_lapl_block(idx) * ao_axis_block(idx) + &
                           2.*(ao_oneD_grad_block_x(idx) * ao_axis_grad_block_x(idx) + &
                               ao_oneD_grad_block_y(idx) * ao_axis_grad_block_y(idx) + &
                               ao_oneD_grad_block_z(idx) * ao_axis_grad_block_z(idx) )
  enddo

END_PROVIDER


BEGIN_PROVIDER [ logical, reduce_primitives ]
 implicit none
 BEGIN_DOC
! If True, remove cusp AOs from AO basis
 END_DOC
 reduce_primitives = do_nucl_fitcusp
! if (calc_density_fit) then
!   reduce_primitives = .False.
! endif
END_PROVIDER


BEGIN_PROVIDER [ logical, primitives_reduced ]
  implicit none
  BEGIN_DOC  
! Tells if the number of primitives has been reduced due to the nucl_fitcusp
  END_DOC
  integer, save :: first_pass = 0

  if (reduce_primitives.and.(first_pass == 0)) then
    first_pass = 1
    integer :: i,j,k,l
    integer :: prim_num_old, prim_num
    prim_num_old = sum(ao_prim_num)
    do k=1,nucl_num
      point(1) = nucl_coord(k,1)
      point(2) = nucl_coord(k,2)
      point(3) = nucl_coord(k,3)+ nucl_fitcusp_radius(k)
      TOUCH point
      PROVIDE  ao_oned_prim_p
      PROVIDE  ao_prim_num_max
      PROVIDE  ao_oned_p
      PROVIDE  ao_power
      PROVIDE  ao_coef
      PROVIDE  ao_nucl
      PROVIDE  mo_fitcusp_normalization_before
      do i=1,ao_num
        if (ao_oned_p(i) /= 0.) then
          l=ao_power(i,1)+ao_power(i,2)+ao_power(i,3)
          if ( (l == 0).and.(ao_nucl(i) == k) ) then
            do j=1,ao_prim_num_max
              if (abs(ao_coef(i,j)*ao_oneD_prim_p(i,j)/ao_oneD_p(i)) < 1.e-4) then
                ao_coef(i,j) = 0.
              endif
            enddo
          endif
        endif
      enddo
    enddo

    do i=1,ao_num
      j=1
      do while (j<ao_prim_num(i))
        if (ao_coef(i,j) == 0.) then
          do k=j,ao_prim_num(i)-1
            ao_coef(i,k) = ao_coef(i,k+1)
            ao_expo(i,k) = ao_expo(i,k+1)
          enddo
          ao_expo(i,ao_prim_num(i)) = 0.
          ao_coef(i,ao_prim_num(i)) = 0.
          ao_prim_num(i) = ao_prim_num(i)-1
        else
          j = j+1
        endif
      enddo
    enddo
    TOUCH ao_expo ao_coef ao_prim_num
    prim_num = sum(ao_prim_num)
    character*(80) :: message
    write(message,'(A8,I4,A)') 'Removed ',prim_num_old-prim_num,' primitive gaussians'
    call info(irp_here,message)
    call iinfo(irp_here,"sum_ao_prim_num",sum(ao_prim_num))
    primitives_reduced = .True.
    FREE ao_oneD_p ao_oneD_prim_p
  else
    primitives_reduced = .False.
  endif

END_PROVIDER



