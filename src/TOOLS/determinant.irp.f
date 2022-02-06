subroutine determinant(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l

 integer :: i,j
 select case (na)
  case default
!DIR$ forceinline
    call determinant_general(a,LDA,na,det_l)
  case (5)
!DIR$ forceinline
    call determinant5(a,LDA,na,det_l)
  case (4)
!DIR$ forceinline
    call determinant4(a,LDA,na,det_l)
  case (3)
!DIR$ forceinline
    call determinant3(a,LDA,na,det_l)
  case (2)
!DIR$ forceinline
    call determinant2(a,LDA,na,det_l)
  case (1)
!DIR$ forceinline
    call determinant1(a,LDA,na,det_l)
  case (0)
    det_l=1.d0
 end select
end

subroutine determinant_general(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l

 double precision :: work(LDA*max(na,64))
!DIR$ ATTRIBUTES ALIGN: $IRP_ALIGN :: work
 integer          :: inf
 integer          :: ipiv(LDA)
!DIR$ ATTRIBUTES ALIGN: $IRP_ALIGN :: ipiv
  integer          :: lwork
  double precision :: f
  integer          :: i, j
  call dgetrf(na, na, a, LDA, ipiv, inf )
  det_l = 1.d0
  j=0
  !DIR$ VECTOR ALIGNED
  do i=1,na
   j = j+min(abs(ipiv(i)-i),1)
   det_l = det_l*a(i,i)
  enddo
  if (iand(j,1) /= 0)  then
    det_l = -det_l
  endif
end

subroutine sdeterminant(a,LDA,na,det_l)
 implicit none
 real             :: a (LDA,na)
 integer          :: LDA
 integer          :: na
 real             :: det_l

 real             :: work(LDA*max(na,64))
!DIR$ ATTRIBUTES ALIGN: $IRP_ALIGN :: work
 integer          :: inf
 integer          :: ipiv(LDA)
!DIR$ ATTRIBUTES ALIGN: $IRP_ALIGN :: ipiv
  integer          :: lwork
  real             :: f
  integer          :: i, j
  call sgetrf(na, na, a, LDA, ipiv, inf )
  det_l = 1.d0
  j=0
  !DIR$ VECTOR ALIGNED
  do i=1,na
   if (ipiv(i) /= i) then
     j = j+1
   endif
   det_l = det_l*a(i,i)
  enddo
  if (iand(j,1) /= 0)  then
    det_l = -det_l
  endif
end

subroutine determinant1(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l

 det_l = a(1,1)
end

subroutine determinant2(a,LDA,na,det_l)
 implicit none
 double precision :: a (LDA,na)
 integer          :: LDA
 integer          :: na
 double precision :: det_l
 double precision :: b(2,2)

 double precision :: f
 det_l = a(1,1)*a(2,2) - a(1,2)*a(2,1)
end

subroutine determinant3(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l
 double precision :: b(4,3)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: b
 integer :: i
 double precision :: f
 det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
        +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

end

subroutine determinant4(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l
 double precision :: b(4,4)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: b
 integer :: i,j
 double precision :: f
 det_l =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                 -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                 +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) &
         -a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                 -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                 +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))) &
         +a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                 -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                 +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) &
         -a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))  &
                 -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))  &
                 +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

end

subroutine determinant5(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l
 double precision :: b(5,5)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: b
 integer :: i,j
 double precision :: f
 det_l = a(1,1)*(a(2,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
 a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))- &
 a(2,3)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)- &
 a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(2,4)*(a(3,2)*( &
 a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+ &
 a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,5)*(a(3,2)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)* &
 a(5,3)-a(4,3)*a(5,2))))-a(1,2)*(a(2,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3)))-a(2,3)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+ &
 a(2,4)*(a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)- &
 a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(2,5)*(a(3,1)*( &
 a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+ &
 a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))+a(1,3)*(a(2,1)*(a(3,2)*(a(4,4)* &
 a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*( &
 a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1)))+a(2,4)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))- &
 a(2,5)*(a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))-a(1,4)*(a(2,1)*( &
 a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,3)* &
 a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1)))-a(2,5)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))+ &
 a(1,5)*(a(2,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*( &
 a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)* &
 a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*( &
 a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,4)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1))))


end

