
! Functions
! =========

subroutine zmq_register_worker(msg)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Register a new worker to the forwarder
  END_DOC
  integer(ZMQ_PTR)               :: msg
  integer                        :: i,rc

  rc = f77_zmq_msg_init_size(msg,8)
  rc = f77_zmq_msg_copy_to_data(msg, 'register',8)
  rc = f77_zmq_msg_send(msg,zmq_to_dataserver_socket,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)


  character*(64) :: buffer
  integer :: size

  size = len(trim(hostname))
  buffer = trim(hostname)
  rc = f77_zmq_msg_init_size(msg,size)
  rc = f77_zmq_msg_copy_to_data(msg, buffer, size)
  rc = f77_zmq_msg_send(msg,zmq_to_dataserver_socket,ZMQ_SNDMORE)
  if (rc == -1) then
     call abrt(irp_here, 'Unable to send register message (1)')
  endif
  rc = f77_zmq_msg_close(msg)

  call worker_log(irp_here, 'Registering')

  rc = f77_zmq_msg_init_size(msg,len_current_PID)
  rc = f77_zmq_msg_copy_to_data(msg, current_PID, len_current_PID)
  rc = f77_zmq_msg_send(msg,zmq_to_dataserver_socket,0)
  if (rc == -1) then
     call abrt(irp_here, 'Unable to send register message (2)')
  endif
  rc = f77_zmq_msg_close(msg)

  rc = f77_zmq_recv(zmq_to_dataserver_socket, buffer, 32, 0)
  if (buffer(1:2)/='OK') then
    call abrt(irp_here, 'Register failed '//trim(http_server))
  endif
  call worker_log(irp_here, 'Registered')

end


subroutine zmq_unregister_worker(msg)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Unregister a new worker to the forwarder
  END_DOC
  integer(ZMQ_PTR)               :: msg
  integer                        :: i,rc

  call worker_log(irp_here, 'Unregistering')
  rc = f77_zmq_msg_init_size(msg,10)
  rc = f77_zmq_msg_copy_to_data(msg, 'unregister',10)
  rc = f77_zmq_msg_send(msg,zmq_to_dataserver_socket,ZMQ_SNDMORE)
  if (rc == -1) then
     call abrt(irp_here, 'Unable to send unregister message (1)')
  endif
  rc = f77_zmq_msg_close(msg)

  character*(64)                 :: buffer
  integer                        :: size

  size = len(trim(hostname))
  buffer = trim(hostname)
  rc = f77_zmq_msg_init_size(msg,size)
  rc = f77_zmq_msg_copy_to_data(msg, buffer, size)
  rc = f77_zmq_msg_send(msg,zmq_to_dataserver_socket,ZMQ_SNDMORE)
  if (rc == -1) then
     call abrt(irp_here, 'Unable to send unregister message (2)')
  endif
  rc = f77_zmq_msg_close(msg)

  rc = f77_zmq_msg_init_size(msg,len_current_PID)
  rc = f77_zmq_msg_copy_to_data(msg, current_PID, len_current_PID)
  rc = f77_zmq_msg_send(msg,zmq_to_dataserver_socket,0)
  if (rc == -1) then
     call abrt(irp_here, 'Unable to send unregister message (3)')
  endif
  rc = f77_zmq_msg_close(msg)

  ! Timeout 15 seconds
  rc = -1
  do i=1,20
    rc = f77_zmq_recv(zmq_to_dataserver_socket, buffer, 32, ZMQ_NOBLOCK)
    if (rc == 2) then
       call worker_log(irp_here, 'Unregistered')
       return
    endif
    call worker_log(irp_here, 'Unregister failed. Retrying')
    call sleep(5)
  enddo
  call abrt(irp_here, 'Unregister failed')


end



subroutine zmq_send_header(msg,header,block_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Receive the header of the multi-part message
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: msg
  character*(*), intent(in)      :: header
  integer, intent(in)            :: block_id
  integer                        :: rc, size
  character*(16)                 :: pid_str

  size = len(trim(header))
  rc = f77_zmq_msg_init_size(msg,size)
  rc = f77_zmq_msg_copy_to_data(msg, header,size)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)

  character*(64) :: buffer

  buffer = trim(hostname)
  size = len(trim(hostname))
  rc = f77_zmq_msg_init_size(msg,size)
  rc = f77_zmq_msg_copy_to_data(msg, buffer, size)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)

  call worker_log(irp_here, header)
  rc = f77_zmq_msg_init_size(msg,len_current_PID)
  rc = f77_zmq_msg_copy_to_data(msg, current_PID, len_current_PID)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)

   
  write(buffer,'(I8)') block_id
  buffer = adjustl(trim(buffer))
  size = len(trim(buffer))
  rc = f77_zmq_msg_init_size(msg,size)
  rc = f77_zmq_msg_copy_to_data(msg, buffer, size)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)

end



subroutine zmq_send_scalar_prop(msg,weight,value)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Send a double precision average over the trajectory
  END_DOC
  integer(ZMQ_PTR)               :: msg
  double precision               :: weight, value
  integer                        :: rc,size
  character*(32)                 :: buffer

  write(buffer,'(E32.16)') weight
  buffer = adjustl(trim(buffer))
  size = len(trim(buffer))
  rc = f77_zmq_msg_init_size(msg,len(trim(buffer)))
  rc = f77_zmq_msg_copy_to_data(msg, buffer,size)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)

  write(buffer,'(E32.16)') value
  buffer = adjustl(trim(buffer))
  size = len(trim(buffer))
  rc = f77_zmq_msg_init_size(msg,size)
  rc = f77_zmq_msg_copy_to_data(msg, buffer,size)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,0)
  rc = f77_zmq_msg_close(msg)

  call worker_log(irp_here,'')
end


subroutine zmq_send_array_prop(msg,weight,value,isize)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Send a double precision average over the trajectory
  END_DOC
  integer(ZMQ_PTR)               :: msg
  integer                        :: isize
  double precision               :: weight, value(isize)
  integer                        :: rc,i,l, sze
  character*(32)                 :: buffer

  write(buffer,'(I8)') isize
  buffer = adjustl(trim(buffer))
  l  = len(trim(buffer))
  rc = f77_zmq_msg_init_size(msg,l)
  rc = f77_zmq_msg_copy_to_data(msg, buffer,l)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)
  sze = l

  write(buffer,'(E32.16)') weight
  buffer = adjustl(trim(buffer))
  l  = len(trim(buffer))
  rc = f77_zmq_msg_init_size(msg,l)
  rc = f77_zmq_msg_copy_to_data(msg, buffer,l)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)
  sze += l

  do i=1,isize
    write(buffer,'(E32.16)') value(i)
    buffer = adjustl(trim(buffer))
    l = len(trim(buffer))
    sze += l
    rc = f77_zmq_msg_init_size(msg,l)
    rc = f77_zmq_msg_copy_to_data(msg, buffer,l)
    if (i < isize) then
      rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
    else
      rc = f77_zmq_msg_send(msg,zmq_socket_push,0)
    endif
    rc = f77_zmq_msg_close(msg)
  enddo

  call worker_log(irp_here,'')
end



subroutine zmq_send_info(msg,message)
  implicit none
  BEGIN_DOC
  ! Send an info message to the forwarder
  END_DOC
  integer(ZMQ_PTR)               :: msg
  character*(64)                 :: message
  integer                        :: rc
  integer                        :: isize

  isize = len(trim(message))
  rc = f77_zmq_msg_init_size(msg,isize)
  rc = f77_zmq_msg_copy_to_data(msg, trim(message),isize)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,0)
  rc = f77_zmq_msg_close(msg)
  call worker_log(irp_here, message)

end


subroutine zmq_send_int(msg,value,isize)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Send an integer array of size n to the forwarder
  END_DOC
  integer(ZMQ_PTR)               :: msg
  integer                        :: isize
  integer                        :: value(isize)
  integer                        :: rc,i,l, sze
  character*(32)                 :: buffer

  write(buffer,'(I8)') isize
  buffer = adjustl(trim(buffer))
  l  = len(trim(buffer))
  rc = f77_zmq_msg_init_size(msg,l)
  rc = f77_zmq_msg_copy_to_data(msg, buffer,l)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)
  sze = l

  do i=1,isize
    write(buffer,'(I16)') value(i)
    buffer = adjustl(trim(buffer))
    l = len(trim(buffer))
    sze += l
    rc = f77_zmq_msg_init_size(msg,l)
    rc = f77_zmq_msg_copy_to_data(msg, buffer,l)
    if (i < isize) then
      rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
    else
      rc = f77_zmq_msg_send(msg,zmq_socket_push,0)
    endif
    rc = f77_zmq_msg_close(msg)
  enddo

  call worker_log(irp_here,'')
end


subroutine zmq_send_real(msg,value,isize)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Send a real array of size n to the forwarder
  END_DOC
  integer(ZMQ_PTR)               :: msg
  integer                        :: isize 
  real                           :: value(isize)
  integer                        :: rc,i,l, sze
  character*(32)                 :: buffer

  write(buffer,'(I8)') isize
  buffer = adjustl(trim(buffer))
  l  = len(trim(buffer))
  rc = f77_zmq_msg_init_size(msg,l)
  rc = f77_zmq_msg_copy_to_data(msg, buffer,l)
  rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
  rc = f77_zmq_msg_close(msg)
  sze = l

  do i=1,isize
    write(buffer,'(E32.16)') value(i)
    buffer = adjustl(trim(buffer))
    l = len(trim(buffer))
    sze += l
    rc = f77_zmq_msg_init_size(msg,l)
    rc = f77_zmq_msg_copy_to_data(msg, buffer,l)
    if (i < isize) then
      rc = f77_zmq_msg_send(msg,zmq_socket_push,ZMQ_SNDMORE)
    else
      rc = f77_zmq_msg_send(msg,zmq_socket_push,0)
    endif
    rc = f77_zmq_msg_close(msg)
  enddo

  call worker_log(irp_here,'')
end



subroutine get_running(do_run)
  use f77_zmq
  implicit none
  include '../types.F'
  BEGIN_DOC
! Fetches the 'do_run' information
  END_DOC
  integer                        :: do_run
  integer                        :: rc, timeout
  character*(16)                 :: buffer
  integer(ZMQ_PTR), save         :: pollitem = 0_ZMQ_PTR

  if (.not.is_worker) then
    do_run = t_Running
    return
  else

    timeout = 600 ! seconds
   
    ! Polling items
    ! -------------
   
    if (pollitem == 0_ZMQ_PTR) then
      pollitem = f77_zmq_pollitem_new(1)
      rc = f77_zmq_pollitem_set_socket(pollitem,1,zmq_socket_running)
      rc = f77_zmq_pollitem_set_events(pollitem,1,ZMQ_POLLIN)
    endif
   
    ! Check for disconnected forwarder after timeout 
    ! ----------------------------------------------
   
    buffer = 'Stopped'
    do while (timeout > 0)
     rc = f77_zmq_poll(pollitem, 1, 100_8)
     if (iand(f77_zmq_pollitem_revents(pollitem,1), ZMQ_POLLIN) /= 0) then
       exit
     endif
     timeout = timeout-1
     call sleep(1)
    enddo
    
    ! Empty the queue to get only the last value
    ! ------------------------------------------
   
    do
      rc = f77_zmq_poll(pollitem, 1, 0)
      if (iand(f77_zmq_pollitem_revents(pollitem,1), ZMQ_POLLIN) == 0) then
         exit
      endif
      rc = f77_zmq_recv(zmq_socket_running, buffer, 16, 0)
    enddo
   
    if (buffer == 'Running') then
      do_run = t_Running
    else if (buffer == 'Queued') then
      do_run = t_Running
    else
      do_run = t_Stopped
    endif
!    call worker_log(irp_here,buffer)

  endif
end
