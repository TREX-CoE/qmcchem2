use f77_zmq



subroutine get_elec_coord_full(elec_coord_full_out,lda)
  implicit none
  integer, intent(in) :: lda
  real, intent(out) :: elec_coord_full_out(lda,3,walk_num)
  integer                        :: rc
  integer                        :: i,j,k,l
  character*(16)                 :: n

  rc = f77_zmq_send(zmq_to_dataserver_socket, 'get_walkers', 11, ZMQ_SNDMORE)
  write(n,*) walk_num
  n = trim(n)
  rc = f77_zmq_send(zmq_to_dataserver_socket, n, len(n), 0)
  call worker_log(irp_here, 'Requesting walkers')

  integer(ZMQ_PTR)               :: msg
  character*(32)                 :: buffer
  integer                        :: sze
  integer                        :: block_time_int

  block_time_int = int(block_time)
  msg = f77_zmq_msg_new()
  sze = 0
  do k=1,walk_num
    do j=1,3
      do i=1,elec_num+1
        rc = f77_zmq_msg_init(msg)
        rc = -1
        do l=1,2*block_time_int
          rc = f77_zmq_msg_recv(msg,zmq_to_dataserver_socket,ZMQ_NOBLOCK)
          if (rc > 0) then
            exit
          endif
          call sleep(1)
        enddo
        if (l>=2*block_time_int) then
          call abrt(irp_here, 'Unable to get walkers')
        endif
        sze += rc
        buffer = ''
        rc = f77_zmq_msg_copy_from_data(msg, buffer)
        rc = f77_zmq_msg_close(msg)
        buffer = trim(adjustl(buffer))
        read(buffer, '(F20.14)') elec_coord_full_out(i,j,k) 
      enddo
    enddo
  enddo
  call worker_log(irp_here, 'Walkers received')
  rc = f77_zmq_msg_destroy(msg)
end





