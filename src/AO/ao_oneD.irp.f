BEGIN_PROVIDER [ logical, compute_all_aos ]
  implicit none
  BEGIN_DOC  
! If true, compute all AOs
  END_DOC

  compute_all_aos = .not. do_nucl_fitcusp
END_PROVIDER


BEGIN_PROVIDER [ integer, ao_elec ]
  implicit none
  BEGIN_DOC  
! Current electron in AO evaluation
  END_DOC
  ao_elec = 1
END_PROVIDER


 BEGIN_PROVIDER [ real, ao_radius, (ao_num_8) ] 
&BEGIN_PROVIDER [ real, ao_radius_sorted, (ao_num_8) ] 
&BEGIN_PROVIDER [ real, nucl_radius, (nucl_num_8) ] 
&BEGIN_PROVIDER [ integer, ao_radius_order, (ao_num_8) ] 
  implicit none
  BEGIN_DOC  
! Radius of the AOs
!
! Radius of the nuclei
!
! Indices of the AOs per nucleus ordered by ao_radius in descending order
  END_DOC
  integer                        :: i,k, inucl
  real                           :: to_sort(ao_num)
  !DIR$ VECTOR ALIGNED
  nucl_radius = 0.
  !DIR$ VECTOR ALIGNED
  do i=1,ao_num
    ao_radius_order(i) = i
  enddo
  do i=1,ao_num
    ao_radius(i) = 0.
    do k=1,ao_prim_num(i)
      ao_radius(i) = max(20./ao_expo_transp(k,i), ao_radius(i))
    enddo
    to_sort(i) = 1./ao_radius(i)
    inucl = ao_nucl(i)
    nucl_radius(inucl) = max(nucl_radius(inucl),ao_radius(i))
  enddo
  !DIR$ VECTOR ALIGNED
  nucl_radius = sqrt(nucl_radius)
  integer                        :: istart, iend
  do inucl=1,nucl_num
    istart = ao_nucl_idx(1,inucl)
    iend   = ao_nucl_idx(2,inucl)
    call isort(to_sort(istart),ao_radius_order(istart),(iend-istart+1))
  enddo
  ao_radius_sorted = ao_radius
  call set_order(ao_radius_sorted,ao_radius_order,ao_num)
     
END_PROVIDER

 BEGIN_PROVIDER [ real, ao_oneD_block, (ao_num_8) ]
&BEGIN_PROVIDER [ integer, ao_oneD_prim_non_zero_idx, ((-simd_sp+1):ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_oneD_lapl_block, (ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_oneD_grad_block_z, (ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_oneD_grad_block_y, (ao_num_8) ]
&BEGIN_PROVIDER [ real, ao_oneD_grad_block_x, (ao_num_8) ]
&BEGIN_PROVIDER [ integer, ao_block_num ]
&BEGIN_PROVIDER [ integer, ao_block_num_8 ]
  implicit none
  BEGIN_DOC
! Gradients and Laplacians of the primitive AOs
  END_DOC
  
  PROVIDE nucl_elec_dist
  PROVIDE ao_power
  PROVIDE ao_expo_unique
  PROVIDE ao_expo_transp
  PROVIDE ao_coef_transp
  PROVIDE compute_all_aos
  PROVIDE ao_prim_num
  PROVIDE ao_expo
  PROVIDE nucl_fitcusp_radius
  
  integer, save                  :: ifirst = 0
  if (ifirst == 0) then
    ifirst = 1
    !DIR$ VECTOR ALIGNED
    ao_oneD_prim_non_zero_idx = 0.
    ao_oneD_block = 0.
    ao_oneD_lapl_block = 0.
    ao_oneD_grad_block_x = 0.
    ao_oneD_grad_block_y = 0.
    ao_oneD_grad_block_z = 0.
  endif
  
  real                           :: ao_oneD_prim_block
  real                           :: buffer2
  
  real                           :: b1, b2, b3, b4
  integer                        :: i, j, k, idx, l, m
  real                           :: rtemp, r, r2
  
  integer                        :: inucl, istart, iend, idx0, ii
  real                           :: buffer(-simd_sp+1:max(ao_nucl_idx(2,nucl_num), ao_expo_unique_nucl_idx(2,nucl_num)) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: buffer
  
  j = ao_elec
  
  real                           :: d
  real, parameter                :: rdm = 20.
  idx = 1
  !DIR$ VECTOR ALIGNED
  buffer = 0.
  do inucl=1,nucl_num
    r = nucl_elec_dist(inucl,j)
    if (r > nucl_radius(inucl)) then
      cycle
    else
      r2 = r*r
      istart = ao_expo_unique_nucl_idx(1,inucl)
      if (ao_expo_unique(istart)*r2 > rdm) then
        cycle
      endif
      real                           :: dcomp
      dcomp = rdm/r2
      b1 = -nucl_elec_dist_vec(1,inucl,j)-nucl_elec_dist_vec(1,inucl,j)
      b2 = -nucl_elec_dist_vec(2,inucl,j)-nucl_elec_dist_vec(2,inucl,j)
      b3 = -nucl_elec_dist_vec(3,inucl,j)-nucl_elec_dist_vec(3,inucl,j)
      b4 =  4.*r2
      iend = ao_expo_unique_nucl_idx(2,inucl)
      
      integer                        :: imax, istep
      imax=istart+1
      !DIR$ LOOP COUNT(20)
      do while (imax<=iend .and. ao_expo_unique(imax) <= dcomp)
        imax += 1
      enddo
      imax = imax-1
      
      buffer(istart) = exp(-r2*ao_expo_unique(istart))
      if (imax - istart == 0) then
        continue
      else if (imax - istart == 1) then
        buffer(istart+1) = exp(-r2*ao_expo_unique(istart+1))
      else if (imax - istart == 2) then
        buffer(istart+1) = exp(-r2*ao_expo_unique(istart+1))
        buffer(istart+2) = exp(-r2*ao_expo_unique(istart+2))
      else
        !DIR$ VECTOR UNALIGNED
        !DIR$ LOOP COUNT MIN(3)
        do i=istart+1,imax
          buffer(i) = exp(-r2*ao_expo_unique(i))
        enddo
      endif
      istart = ao_nucl_idx(1,inucl)
      iend   = ao_nucl_idx(2,inucl)
      idx0 = idx
      
      if ( r > nucl_fitcusp_radius(inucl) ) then
        
        !DIR$ VECTOR UNALIGNED
        do ii=istart,iend
          if ( r2 > ao_radius_sorted(ii)) then
            exit
          endif
          ao_oneD_prim_non_zero_idx(idx) = ao_radius_order(ii)
          idx += 1
        enddo
        
      else if (.not.compute_all_aos) then
        
        !DIR$ VECTOR UNALIGNED
        do ii=istart,iend
          if ( r2 > ao_radius_sorted(ii)) then
            exit
          endif
          i = ao_radius_order(ii)
          if ( .not.ao_power_is_zero(i) ) then
            ao_oneD_prim_non_zero_idx(idx) = i
            idx += 1
          endif
        enddo
        
      else
        
        !DIR$ VECTOR UNALIGNED
        do ii=istart,iend
          if ( r2 > ao_radius_sorted(ii)) then
            exit
          endif
          ao_oneD_prim_non_zero_idx(idx) = ao_radius_order(ii)
          idx += 1
        enddo
        
      endif
      
      
      iend = idx-1
      
      real                           :: f, g, h
      integer                        :: lmax
      i = ao_oneD_prim_non_zero_idx(idx0)
      if (i==0) then
        cycle
      endif
      lmax = ao_prim_num(i)
      g = buffer(ao_expo_unique_idx_1(i))
      f = ao_coef(i,1)*g*ao_expo(i,1)
      ao_oneD_block(idx0) = ao_coef(i,1)*g
      ao_oneD_lapl_block(idx0) = f*(b4 * ao_expo(i,1)-6.)
      ao_oneD_grad_block_x(idx0) = b1*f
      ao_oneD_grad_block_y(idx0) = b2*f
      ao_oneD_grad_block_z(idx0) = b3*f
      
      select case (iend-idx0)
        case (0)
          continue
        case (1)
          i = ao_oneD_prim_non_zero_idx(idx0+1)
          g = buffer(ao_expo_unique_idx_1(i))
          f = ao_coef(i,1)*g*ao_expo(i,1)
          ao_oneD_block(idx0+1) = ao_coef(i,1)*g
          ao_oneD_lapl_block(idx0+1) = f*(b4 * ao_expo(i,1)-6.)
          ao_oneD_grad_block_x(idx0+1) = b1*f
          ao_oneD_grad_block_y(idx0+1) = b2*f
          ao_oneD_grad_block_z(idx0+1) = b3*f
          lmax = max(ao_prim_num(i),lmax)
        case (2)
          i = ao_oneD_prim_non_zero_idx(idx0+1)
          g = buffer(ao_expo_unique_idx_1(i))
          f = ao_coef(i,1)*g*ao_expo(i,1)
          ao_oneD_block(idx0+1) = ao_coef(i,1)*g
          ao_oneD_lapl_block(idx0+1) = f*(b4 * ao_expo(i,1)-6.)
          ao_oneD_grad_block_x(idx0+1) = b1*f
          ao_oneD_grad_block_y(idx0+1) = b2*f
          ao_oneD_grad_block_z(idx0+1) = b3*f
          lmax = max(ao_prim_num(i),lmax)
          
          i = ao_oneD_prim_non_zero_idx(idx0+2)
          g = buffer(ao_expo_unique_idx_1(i))
          f = ao_coef(i,1)*g*ao_expo(i,1)
          ao_oneD_block(idx0+2) = ao_coef(i,1)*g
          ao_oneD_lapl_block(idx0+2) = f*(b4 * ao_expo(i,1)-6.)
          ao_oneD_grad_block_x(idx0+2) = b1*f
          ao_oneD_grad_block_y(idx0+2) = b2*f
          ao_oneD_grad_block_z(idx0+2) = b3*f
          lmax = max(ao_prim_num(i),lmax)
        case (3)
          i = ao_oneD_prim_non_zero_idx(idx0+1)
          g = buffer(ao_expo_unique_idx_1(i))
          f = ao_coef(i,1)*g*ao_expo(i,1)
          ao_oneD_block(idx0+1) = ao_coef(i,1)*g
          ao_oneD_lapl_block(idx0+1) = f*(b4 * ao_expo(i,1)-6.)
          ao_oneD_grad_block_x(idx0+1) = b1*f
          ao_oneD_grad_block_y(idx0+1) = b2*f
          ao_oneD_grad_block_z(idx0+1) = b3*f
          lmax = max(ao_prim_num(i),lmax)
          
          i = ao_oneD_prim_non_zero_idx(idx0+2)
          g = buffer(ao_expo_unique_idx_1(i))
          f = ao_coef(i,1)*g*ao_expo(i,1)
          ao_oneD_block(idx0+2) = ao_coef(i,1)*g
          ao_oneD_lapl_block(idx0+2) = f*(b4 * ao_expo(i,1)-6.)
          ao_oneD_grad_block_x(idx0+2) = b1*f
          ao_oneD_grad_block_y(idx0+2) = b2*f
          ao_oneD_grad_block_z(idx0+2) = b3*f
          lmax = max(ao_prim_num(i),lmax)
          
          i = ao_oneD_prim_non_zero_idx(idx0+3)
          g = buffer(ao_expo_unique_idx_1(i))
          f = ao_coef(i,1)*g*ao_expo(i,1)
          ao_oneD_block(idx0+3) = ao_coef(i,1)*g
          ao_oneD_lapl_block(idx0+3) = f*(b4 * ao_expo(i,1)-6.)
          ao_oneD_grad_block_x(idx0+3) = b1*f
          ao_oneD_grad_block_y(idx0+3) = b2*f
          ao_oneD_grad_block_z(idx0+3) = b3*f
          lmax = max(ao_prim_num(i),lmax)
          
          case default
          !DIR$ VECTOR UNALIGNED
          !DIR$ LOOP COUNT MIN(4)
          do idx=idx0+1,iend
            i = ao_oneD_prim_non_zero_idx(idx)
            g = buffer(ao_expo_unique_idx_1(i))
            f = ao_coef(i,1)*g*ao_expo(i,1)
            ao_oneD_block(idx) = ao_coef(i,1)*g
            ao_oneD_lapl_block(idx) = f*(b4 * ao_expo(i,1)-6.)
            ao_oneD_grad_block_x(idx) = b1*f
            ao_oneD_grad_block_y(idx) = b2*f
            ao_oneD_grad_block_z(idx) = b3*f
            lmax = max(ao_prim_num(i),lmax)
          enddo

      end select
      
      if (lmax == 1) then
        idx = iend+1
        cycle
      endif
      
      do idx=idx0,iend
        i = ao_oneD_prim_non_zero_idx(idx)
        if (ao_prim_num(i) == 1) then
          cycle
        endif
        l = ao_prim_num(i)
        
        ao_oneD_block(idx) = 0.
        ao_oneD_lapl_block(idx) = 0.
        f = 0.
        g = 0.
        h = 0.
        !DIR$ VECTOR ALIGNED
        !DIR$ LOOP COUNT (8)
        do k=1,l
          ao_oneD_prim_block = buffer(ao_expo_unique_idx(k,i))*ao_coef_transp(k,i)
          g = g + ao_oneD_prim_block
          buffer2 = ao_expo_transp(k,i)*ao_oneD_prim_block
          f = f+buffer2
          h = h + buffer2 * (b4*ao_expo_transp(k,i)-6.)
        enddo
        ao_oneD_block(idx) = g
        ao_oneD_lapl_block(idx) = h
        ao_oneD_grad_block_x(idx) = b1*f
        ao_oneD_grad_block_y(idx) = b2*f
        ao_oneD_grad_block_z(idx) = b3*f
      enddo
      idx = iend+1
    endif
  enddo
  
  ao_oneD_prim_non_zero_idx(0) = idx-1
  if (ao_oneD_prim_non_zero_idx(0) > ao_block_num) then
    ao_block_num = ao_oneD_prim_non_zero_idx(0)
    integer                        :: mod_align
    ao_block_num_8 = mod_align(ao_block_num)
  endif
  
  if (ao_oneD_prim_non_zero_idx(0) == 0) then
    ao_oneD_prim_non_zero_idx(0) = 1
    ao_oneD_prim_non_zero_idx(1) = 1
    ao_oneD_block(idx) = 0.
    ao_oneD_lapl_block(idx) = 0.
    ao_oneD_grad_block_x(idx) = 0.
    ao_oneD_grad_block_y(idx) = 0.
    ao_oneD_grad_block_z(idx) = 0.
  endif
  
END_PROVIDER


