let _list    = ref []
let _running = ref false
let _threads = ref []

(** Kill the current process and all children *)
let kill () =
  let kill pid =
    Unix.kill pid Sys.sigkill;
    Printf.printf "Killed %d\n%!" pid
  in
  List.iter kill (!_list);
  exit 1



(** Start watchdog *)
let start () =

  if (!_running) then
    failwith "Watchdog error: Already running"
  else
    begin
      _running := true;

      let pause () =
        Unix.sleep 1
      in

      let pid_is_running pid =
         Sys.file_exists ("/proc/"^(string_of_int pid)^"/stat")
      in

      let f () =
        while (!_running)
        do
          pause () ;

(*DEBUG
 List.iter (fun x -> Printf.printf "%d\n%!" x) (!_list) ;
 *)

          let continue () =
            List.fold_left
              ( fun accu x -> accu && (pid_is_running x))
              true (!_list)
          in
          if ( not (continue ()) ) then
            kill ()
        done
      in
      _threads := ( (Thread.create f) () ) :: (!_threads)
    end


(** Stop watchdog *)
let stop () =
  if (!_running) then
    _running := false
  else
    failwith "Watchdog error: Already stopped"


(** Add a PID to tracking *)
let add pid =
  if (not !_running) then
    start ();
  _list := pid :: (!_list)


(** Remove a PID from tracking *)
let del pid =
  let rec aux accu = function
    | [] -> accu
    | a :: rest ->
        if (a <> pid) then
          aux (a::accu) rest
        else
          aux accu rest
  in
  _list := aux [] (!_list);

  match (!_list) with
  | [] -> if (!_running) then stop ()
  | _ -> ()


(** Fork and exec a new process *)
let fork_exec ~prog ~args () =
  let pid =
    match Unix.fork () with
    | 0 -> Unix.execvp prog args
    | pid -> pid
  in

  let f () =
    add pid;
    let success =
      match (Unix.waitpid [] pid) with
      | pid , Unix.WEXITED n -> true
      | pid , Unix.WSIGNALED n ->
          ( Printf.printf "PID %d killed with signal %d\n%!" pid n;
            false )
      | pid , Unix.WSTOPPED n ->
          ( Printf.printf "PID %d stopped with signal %d\n%!" pid n;
            false )
    in
    del pid ;
    if (not success) then
      kill ()
  in
  _threads := ( (Thread.create f) () ) :: (!_threads);
  pid


(** Wait for threads to finish *)
let join () =
(*  if (!_running) then stop (); *)
  List.iter Thread.join (!_threads);
  assert (not !_running)

