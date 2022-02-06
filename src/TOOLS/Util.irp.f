 BEGIN_PROVIDER [ character*(8), current_PID ]
&BEGIN_PROVIDER [ integer,   len_current_PID ]
 implicit none
 BEGIN_DOC
 ! Process ID
 END_DOC
 integer :: getpid
 write(current_PID,'(I8)') getpid()
 current_PID = adjustl(trim(current_PID))
 len_current_PID = len(trim(current_PID))

END_PROVIDER


 BEGIN_PROVIDER [ integer, simd_sp ]
&BEGIN_PROVIDER [ integer, simd_dp ]
 implicit none
 BEGIN_DOC  
! Number of array elements for vectorization
 END_DOC
 simd_sp = max(1,$IRP_ALIGN / 4)
 simd_dp = max(1,$IRP_ALIGN / 8)
END_PROVIDER

integer function mod_align(n)
 implicit none
 integer, intent(in) :: n

 if (mod(n,simd_sp) /= 0) then
   mod_align = n + simd_sp - mod(n,simd_sp)
 else
   mod_align = n
 endif
end function

real function fact2(n)
 implicit none
 integer :: n
 real :: dblefact_even
 real :: dblefact_odd
 if (mod(n,2) == 0) then
   fact2 = dblefact_even(n)
 else
   fact2 = dblefact_odd(n)
 endif
end


real function dblefact_even(n)
  implicit none
  integer :: n
  real, save :: memo(1:100)
  integer, save :: memomax = 2
  integer :: i
  if (n<=memomax) then
    if (n<2) then
      dblefact_even = 1.
    else
      dblefact_even = memo(n)
    endif
    return
  endif

  memo(2) = 2.
  do i=memomax+2,min(n,100),2
    memo(i) = memo(i-2)* float(i)
  enddo
  memomax = min(n,100)
  dblefact_even = memo(memomax)

  do i=102,n,2
    dblefact_even = dblefact_even*float(i)
  enddo

end


real function dblefact_odd(n)
  implicit none
  integer :: n
  real, save :: memo(1:100)
  integer, save :: memomax = 1
  integer :: i

  if (n<=memomax) then
    if (n<3) then
      dblefact_odd = 1.
    else
      dblefact_odd = memo(n)
    endif
    return
  endif

  memo(1) = 1.
  do i=memomax+2,min(n,99),2
    memo(i) = memo(i-2)* float(i)
  enddo
  memomax = min(n,99)
  dblefact_odd = memo(memomax)

  do i=101,n,2
    dblefact_odd = dblefact_odd*float(i)
  enddo

end



real function fact(n)
  implicit none
  integer :: n
  real , save :: memo(1:100)
  integer, save :: memomax = 1 
  
  if (n<=memomax) then
    if (n<2) then
      fact = 1.
    else
      fact = memo(n)
    endif
    return
  endif

  integer :: i
  memo(1) = 1.
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)*float(i)
  enddo
  memomax = min(n,100)
  fact = memo(memomax)
  do i=101,n
    fact = fact*float(i)
  enddo
end function


real function rintgauss(n)
  implicit none
  include '../constants.F'

  integer :: n

  rintgauss = sqrt(pi)
  if ( n == 0 ) then
    return
  else if ( n == 1 ) then
    rintgauss = 0.
  else if ( mod(n,2) == 1) then
    rintgauss = 0.
  else
    real :: fact2
    rintgauss = rintgauss/(2.**(n/2))
    rintgauss = rintgauss * fact2(n-1)
  endif
end function

real function goverlap(gamA,gamB,nA)
  implicit none

  real    :: gamA, gamB
  integer :: nA(3)

  real :: gamtot
  gamtot = gamA+gamB

  goverlap=1.0

  integer :: l
  real :: rintgauss
  do l=1,3
    goverlap *= rintgauss(nA(l)+nA(l))/ (gamtot**(0.5+float(nA(l))))
  enddo

end function

double precision function binom(n,k)
 implicit none
 integer, intent(in) :: k,n
 real :: fact
 binom=fact(n)/(fact(k)*fact(n-k))
end

!DIR$ attributes forceinline :: transpose
recursive subroutine transpose(A,LDA,B,LDB,d1,d2)
 implicit none
 BEGIN_DOC
! Transpose input matrix A into output matrix B
 END_DOC
 integer, intent(in) :: d1, d2, LDA, LDB
 real, intent(in) :: A(LDA,d2)
 real, intent(out) :: B(LDB,d1)

 integer :: i,j,k, mod_align
 if ( d2 < 32 ) then
   do j=1,d1
     !DIR$ LOOP COUNT (16)
     do i=1,d2
      B(i,j  ) = A(j  ,i)
     enddo
    enddo
   return
 else if (d1 > d2) then
 !DIR$ forceinline 
   k=d1/2
 !DIR$ forceinline recursive
   call transpose(A(1,1),LDA,B(1,1),LDB,k,d2)
 !DIR$ forceinline recursive
   call transpose(A(k+1,1),LDA,B(1,k+1),LDB,d1-k,d2)
   return
 else
 !DIR$ forceinline
   k=d2/2
 !DIR$ forceinline recursive
   call transpose(A(1,k+1),LDA,B(k+1,1),LDB,d1,d2-k)
 !DIR$ forceinline recursive
   call transpose(A(1,1),LDA,B(1,1),LDB,d1,k)
   return
 endif

end

subroutine gaussian_product(a,xa,b,xb,k,p,xp)
 implicit none
! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}

 real, intent(in) :: a,b         ! Exponents
 real, intent(in) :: xa(3),xb(3) ! Centers
 real, intent(out) :: p          ! New exponent
 real, intent(out) :: xp(3)      ! New center
 real, intent(inout) :: k        ! Constant

 real :: p_inv

 ASSERT (a>0.)
 ASSERT (b>0.)

 real :: xab(4), ab
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xab

 p_inv = 1./(a+b)
 p = a+b
 ab = a*b
 xab(1) = xa(1)-xb(1)
 xab(2) = xa(2)-xb(2)
 xab(3) = xa(3)-xb(3)
 xp(1) = (a*xa(1)+b*xb(1))*p_inv
 xp(2) = (a*xa(2)+b*xb(2))*p_inv
 xp(3) = (a*xa(3)+b*xb(3))*p_inv
 k *= exp(-ab*p_inv*(xab(1)*xab(1)+xab(2)*xab(2)+xab(3)*xab(3)))

end subroutine

!DIR$ attributes forceinline :: transpose_to_dp
recursive subroutine transpose_to_dp(A,LDA,B,LDB,d1,d2)
 implicit none
 BEGIN_DOC
! Transpose SP input matrix A into DP output matrix B
 END_DOC
 integer, intent(in) :: d1, d2, LDA, LDB
 real, intent(in) :: A(LDA,d2)
 double precision, intent(out) :: B(LDB,d1)

 integer :: i,j,k, mod_align
 if ( d2 < 32 ) then
   do j=1,d1
     !DIR$ LOOP COUNT (16)
     do i=1,d2
      B(i,j  ) = A(j  ,i)
     enddo
    enddo
   return
 else if (d1 > d2) then
 !DIR$ forceinline 
   k=d1/2
 !DIR$ forceinline recursive
   call transpose_to_dp(A(1,1),LDA,B(1,1),LDB,k,d2)
 !DIR$ forceinline recursive
   call transpose_to_dp(A(k+1,1),LDA,B(1,k+1),LDB,d1-k,d2)
   return
 else
 !DIR$ forceinline
   k=d2/2
 !DIR$ forceinline recursive
   call transpose_to_dp(A(1,k+1),LDA,B(k+1,1),LDB,d1,d2-k)
 !DIR$ forceinline recursive
   call transpose_to_dp(A(1,1),LDA,B(1,1),LDB,d1,k)
   return
 endif

end

