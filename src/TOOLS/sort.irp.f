BEGIN_TEMPLATE
 subroutine insertion_$Xsort (x,iorder,isize)
  implicit none
  $type,intent(inout)    :: x(isize)
  integer,intent(inout)  :: iorder(isize)
  integer,intent(in)     :: isize
  $type                  :: xtmp
  integer                :: i, i0, j, jmax
 
  do i=1,isize
   xtmp = x(i)
   i0 = iorder(i)
   j = i-1
   do j=i-1,1,-1
    if ( x(j) > xtmp ) then 
     x(j+1) = x(j)
     iorder(j+1) = iorder(j)
    else
     exit
    endif
   enddo
   x(j+1) = xtmp
   iorder(j+1) = i0
  enddo
 
 end subroutine insertion_$Xsort

 subroutine heap_$Xsort(x,iorder,isize)
  implicit none
  $type,intent(inout)    :: x(isize)
  integer,intent(inout)  :: iorder(isize)
  integer,intent(in)     :: isize
  
  integer :: i, k, j, l, i0
  $type :: xtemp

 l = isize/2+1
 k = isize
 do while (.True.)
    if (l>1) then
      l=l-1
      xtemp = x(l)
      i0 = iorder(l)
    else
      xtemp = x(k)
      i0 = iorder(k)
      x(k) = x(1)
      iorder(k) = iorder(1)
      k = k-1
      if (k == 1) then
        x(1) = xtemp 
        iorder(1) = i0
        exit
      endif
    endif
    i=l
    j = shiftl(l,1)
    do while (j<k)
     if ( x(j) < x(j+1) ) then
       j=j+1
     endif
     if (xtemp < x(j)) then
       x(i) = x(j)
       iorder(i) = iorder(j)
       i = j
       j = shiftl(j,1)
     else
       j = k+1
     endif
    enddo
    if (j==k) then
     if (xtemp < x(j)) then
       x(i) = x(j)
       iorder(i) = iorder(j)
       i = j
       j = shiftl(j,1)
     else
       j = k+1
     endif
    endif
    x(i) = xtemp 
    iorder(i) = i0
 enddo

 end subroutine heap_$Xsort

 subroutine $Xsort(x,iorder,isize)
  implicit none
  $type,intent(inout)    :: x(isize)
  integer,intent(inout)  :: iorder(isize)
  integer,intent(in)     :: isize
  if (isize < 32) then
    call insertion_$Xsort(x,iorder,isize)
  else
    call heap_$Xsort(x,iorder,isize)
  endif
 end subroutine $Xsort

SUBST [ X, type ]
   ; real ;;
 d ; double precision ;;
 i ; integer ;;
END_TEMPLATE

BEGIN_TEMPLATE
 subroutine $Xset_order(x,iorder,isize)
  implicit none
  integer          :: isize
  $type            :: x(*), xtmp(isize)
  integer          :: iorder(*)
  integer          :: i 
 
  do i=1,isize
   xtmp(i) = x(iorder(i))
  enddo
 
  do i=1,isize
   x(i) = xtmp(i)
  enddo
 end

SUBST [ X, type ]
   ; real ;;
 d ; double precision ;;
 i ; integer ;;
 l ; logical ;;
END_TEMPLATE

