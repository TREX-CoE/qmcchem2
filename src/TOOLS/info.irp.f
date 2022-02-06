BEGIN_TEMPLATE

 subroutine $Xinfo (here,token,value)
  implicit none
  character*(*), intent(in)    :: here
  character*(*), intent(in)    :: token
  $Y, intent(in)         :: value
  if (print_level>1) then
   write(0,*)  trim(here)//':'
   $Z
  endif
 end
 
SUBST [ X, Y, Z ]
  r; real; 
   write(0,*)  ' -> ', trim(token), '=', value;;
  d; double precision; 
   write(0,*)  ' -> ', trim(token), '=', value;;
  i; integer; 
   write(0,*)  ' -> ', trim(token), '=', value;;
  c; character*(*); 
   write(0,*)  ' -> ', trim(token), '=', value;;
  l; logical;
   if (value) then
     write(0,*)  ' -> ', trim(token), '= True'
   else
     write(0,*)  ' -> ', trim(token), '= False'
   endif ;;
END_TEMPLATE

subroutine info(here,message)
  implicit none
  character*(*), intent(in) :: here, message
  if (print_level > 1) then
     write(0,*)  trim(here)//':'
     write(0,*)  ' -> ', trim(message)
  endif
end
