subroutine abrt (here,message)
  implicit none
  character*(*)                  :: here
  character*(*)                  :: message
  write(0,*) '-------------------------'
  write(0,*) 'Error in '//trim(here)//':'
  write(0,*) '-------------------------'
  write(0,*) trim(message)//'.'
  write(0,*) '-------------------------'
  if (is_worker) then
    call worker_log(here,message)
    call sleep(2)
  endif
  call finish()
  stop
end

subroutine finish()
  implicit none
  call ezfio_finish()
end

