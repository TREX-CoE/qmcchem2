subroutine main_qmc
 use f77_zmq
 implicit none

 is_worker = .True.
 SOFT_TOUCH is_worker 
 call start_main_qmc()

end

subroutine start_main_qmc
 use f77_zmq
 implicit none
 integer(ZMQ_PTR)               :: msg
 integer                        :: rc, v
 integer*8                      :: cpu0, count_rate, count_max

 ! Initialization 
 ! --------------

 call system_clock(cpu0, count_rate, count_max)

 msg   = f77_zmq_msg_new()
 call zmq_register_worker(msg)

 ! One equilibration block
 ! -----------------------

 call equilibration

 ! Run the QMC blocks
 ! ------------------

 call run_qmc(cpu0)

 ! Clean exit
 ! ----------

 call zmq_unregister_worker(msg)
 rc = f77_zmq_msg_destroy(msg)
 v = 0
 rc = f77_zmq_setsockopt(zmq_socket_push,ZMQ_LINGER,v,4)
 rc = f77_zmq_setsockopt(zmq_socket_running,ZMQ_LINGER,v,4)
 rc = f77_zmq_setsockopt(zmq_to_dataserver_socket,ZMQ_LINGER,v,4)
 rc = f77_zmq_close(zmq_socket_push)
 rc = f77_zmq_close(zmq_socket_running)
 rc = f77_zmq_close(zmq_to_dataserver_socket)
! rc = f77_zmq_ctx_destroy(zmq_context)
end

subroutine equilibration
 PROVIDE E_loc_block_walk
end





subroutine run_qmc(cpu0)
 use f77_zmq
 implicit none
 include '../types.F'
 
 integer*8                      :: cpu0
 integer                        :: isize, i, j, ierr
 double precision               :: min, max
 real                           :: value
 integer(ZMQ_PTR)               :: msg
 integer                        :: do_run, rc
 integer                        :: block_id
 integer*8                      :: cpu1, count_rate, count_max
 


 PROVIDE elec_num

 msg   = f77_zmq_msg_new()
 call get_running(do_run)

 block_id = 0

 do while (do_run == t_Running)

  block_id += 1
  call accep_reset
  TOUCH elec_coord
  PROVIDE block_weight E_loc_block_walk

  ! Start by sending accept rate
  real, external :: accep_rate
  call zmq_send_header(msg,'accep',block_id)
  value = accep_rate()
  call zmq_send_real(msg, value, 1)

  double precision :: v0, v1, v2, h
  double precision :: d1v,d2v,d1e,d2e


  if (save_data) then
    call zmq_send_header(msg,'elec_coord',block_id)
    call zmq_send_real(msg,elec_coord_full,size(elec_coord_full))
  endif

BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

derivlist = []

td = """
     do j=0,size($X_block_walk,1)-1,7
      if ($X_block_walk(j+Pos_weight) == 0.d0) then
        $X_block_walk(j+1) = 0.d0
        $X_block_walk(j+2) = 0.d0
        $X_block_walk(j+3) = 0.d0
        $X_block_walk(j+4) = 0.d0
        $X_block_walk(j+5) = 0.d0
        $X_block_walk(j+6) = 0.d0
        $X_block_walk(j+7) = 0.d0
        cycle
      endif
      $X_block_walk(j+Pos_E_loc) =    &
         $X_block_walk(j+Pos_E_loc) / $X_block_walk(j+Pos_weight)
      $X_block_walk(j+Pos_E_loc_2) =    &
         $X_block_walk(j+Pos_E_loc_2) / $X_block_walk(j+Pos_weight)

      $X_block_walk(j+Neg_E_loc) =    &
         $X_block_walk(j+Neg_E_loc) / $X_block_walk(j+Neg_weight)
      $X_block_walk(j+Neg_E_loc_2) =    &
         $X_block_walk(j+Neg_E_loc_2) / $X_block_walk(j+Neg_weight)
      h = $X_block_walk(j+Delta)
      v0 = E_loc_block_walk
      v1 = $X_block_walk(j+Pos_E_loc)
      v2 = $X_block_walk(j+Neg_E_loc)
      d1e = 0.5d0*(v1-v2)/h
      d2e = (v1+v2-v0-v0)/(h*h)

      v0 = dabs(E_loc_2_block_walk - v0*v0)
      v1 = dabs($X_block_walk(j+Pos_E_loc_2) - $X_block_walk(j+Pos_E_loc)**2)
      v2 = dabs($X_block_walk(j+Neg_E_loc_2) - $X_block_walk(j+Neg_E_loc)**2)
      d1v = 0.5d0*(v1-v2)/h
      d2v = (v1+v2-v0-v0)/(h*h)

      $X_block_walk(j+1) = d1e
      $X_block_walk(j+2) = d2e
      $X_block_walk(j+3) = d1v
      $X_block_walk(j+4) = d2v
      $X_block_walk(j+5) = 0.d0
      $X_block_walk(j+6) = 0.d0
      $X_block_walk(j+7) = 0.d0
    enddo
"""

for p in properties:
  t = """
    if (calc_$X) then
    """
  if p[2] == "":
    t += """
      call zmq_send_header(msg,'$X',block_id)
      call zmq_send_scalar_prop(msg,block_weight,$X_block_walk)
      $X_2_block_walk = dabs($X_2_block_walk - $X_block_walk*$X_block_walk)
      call zmq_send_header(msg,'$X_qmcvar',block_id)
      call zmq_send_scalar_prop(msg,block_weight,$X_2_block_walk)
    """
  else:
    if p[1] in derivlist:
      t+= td
    t += """
      isize = size($X_block_walk)
      call zmq_send_header(msg,'$X',block_id)
      call zmq_send_array_prop(msg,block_weight,$X_block_walk,isize)
      """
  t += """
! TODO : Min and Max are commented here
!     call zmq_send_header(msg,'$X_min',block_id)
!     call zmq_send_real(msg,$X_min)
!     call zmq_send_header(msg,'$X_max',block_id)
!     call zmq_send_real(msg,block_weight,$X_max)
   endif
  """
  print (t.replace("$X",p[1]))
END_SHELL



  ! Finish by sending CPU time
  call system_clock(cpu1, count_rate, count_max)
  value = real(cpu1-cpu0)/real(count_rate)
  call zmq_send_header(msg,'cpu',block_id)
  call zmq_send_real(msg, value, 1)

  call get_running(do_run)
  cpu0 = cpu1
 enddo
99 continue
!
end
