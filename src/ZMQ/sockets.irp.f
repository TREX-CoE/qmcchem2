use f77_zmq

! Addresses
! =========

BEGIN_PROVIDER [ character*(48), dataserver_address ]
  implicit none
  BEGIN_DOC  
! Address of the data server
  END_DOC
  dataserver_address = trim(http_server)
  integer :: i
  do i=len(dataserver_address),1,-1
    if ( dataserver_address(i:i) == ':') then
      dataserver_address = trim(dataserver_address(1:i-1))
      exit
    endif
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, zmq_port_start ]
 implicit none
 BEGIN_DOC
 ! Starting port for ZMQ
 END_DOC
 zmq_port_start = -1
 double precision :: qmc_ranf
 integer :: i,l
 character*(8) :: buffer
 l = len(http_server)
 do i=len(http_server),1,-1
   if ( http_server(i:i) == ':') then
     buffer = trim(http_server(i+1:l))
     read(buffer, *) zmq_port_start
     exit
   endif
 enddo
END_PROVIDER

function zmq_port(ishift)
  implicit none
  integer, intent(in)            :: ishift
  character*(8)                  :: zmq_port
  write(zmq_port,'(I8)') zmq_port_start+ishift
  zmq_port = adjustl(trim(zmq_port))
end

! Sockets
! =======

BEGIN_PROVIDER [ integer(ZMQ_PTR), zmq_context ]
 implicit none
 BEGIN_DOC
 ! Context for the ZeroMQ library
 END_DOC
 zmq_context = f77_zmq_ctx_new ()
END_PROVIDER

BEGIN_PROVIDER [ integer(ZMQ_PTR), zmq_to_dataserver_socket ]
 implicit none
 BEGIN_DOC
 ! Socket on which the dataserver replies
 END_DOC
 integer                        :: rc
 zmq_to_dataserver_socket = f77_zmq_socket(zmq_context, ZMQ_REQ)
 rc = f77_zmq_connect(zmq_to_dataserver_socket, trim(http_server))
 if (rc /= 0) then
   call abrt(irp_here, 'Unable to connect zmq_to_dataserver_socket')
 endif
 integer                        :: i
 i=4
 rc = f77_zmq_setsockopt(zmq_to_dataserver_socket, ZMQ_SNDTIMEO, 600000, i)
 if (rc /= 0) then
   call abrt(irp_here, 'Unable to set send timout in zmq_to_dataserver_socket')
 endif
 rc = f77_zmq_setsockopt(zmq_to_dataserver_socket, ZMQ_RCVTIMEO, 600000, i)
 if (rc /= 0) then
   call abrt(irp_here, 'Unable to set recv timout in zmq_to_dataserver_socket')
 endif
 call worker_log(irp_here,'REQ socket : '//trim(http_server))
END_PROVIDER

BEGIN_PROVIDER [ integer(ZMQ_PTR), zmq_socket_running ]
 implicit none
 BEGIN_DOC
 ! Socket on which the dataserver sends the running status
 END_DOC
 integer                        :: rc
 character*(64)                 :: address
 character*(8), external        :: zmq_port
 zmq_socket_running = f77_zmq_socket(zmq_context, ZMQ_SUB)
 address = trim(dataserver_address)//':'//zmq_port(1)
 rc = f77_zmq_connect(zmq_socket_running, trim(address))
 if (rc /= 0) then
   call abrt(irp_here, 'Unable to connect zmq_socket_running')
 endif
 rc = f77_zmq_setsockopt(zmq_socket_running,ZMQ_SUBSCRIBE,'',0)
 call worker_log(irp_here,'Running socket : '//trim(address))
END_PROVIDER

BEGIN_PROVIDER [ integer(ZMQ_PTR), zmq_socket_push ]
 implicit none
 BEGIN_DOC
 ! Socket on which to push the results
 END_DOC
 integer                        :: rc
 character*(64)                 :: address
 character*(8), external        :: zmq_port
 zmq_socket_push = f77_zmq_socket(zmq_context, ZMQ_PUSH)
 address = trim(dataserver_address)//':'//zmq_port(2)
 rc = f77_zmq_setsockopt(zmq_socket_push,ZMQ_LINGER,600000,4)
 rc = f77_zmq_connect(zmq_socket_push, trim(address))
 if (rc /= 0) then
   call abrt(irp_here, 'Unable to connect zmq_socket_push')
 endif
 call worker_log(irp_here,'Push socket : '//trim(address))

END_PROVIDER

BEGIN_PROVIDER [ integer(ZMQ_PTR), zmq_socket_log ]
 implicit none
 BEGIN_DOC
 ! Socket on which to send the log
 END_DOC
 integer                        :: rc
 character*(64)                 :: address
 character*(8), external        :: zmq_port
 zmq_socket_log = f77_zmq_socket(zmq_context, ZMQ_PUB)
 address = trim(dataserver_address)//':'//zmq_port(3)
 rc = f77_zmq_connect(zmq_socket_log, trim(address))

END_PROVIDER

subroutine worker_log(where_, message)
  implicit none
  character*(*), intent(in)      :: where_
  character*(*), intent(in)      :: message
  character*(512)                :: buffer
  integer                        :: rc
  if (is_worker) then
    write(buffer,'(A,X,A)') trim(hostname)//':'//trim(current_pid), &
    trim(where_)//' : '//trim(message)
    rc = f77_zmq_send(zmq_socket_log,buffer,len_trim(buffer),0)
  endif
end
