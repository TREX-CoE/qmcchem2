
let full_run ?(start_dataserver=true) ezfio_filename  =
  (* Identify the job scheduler *)
  let launcher =
    Launcher.find ()
  and scheduler =
    Scheduler.find ()
  in
  Printf.printf "Scheduler : %s\n%!" (Scheduler.to_string scheduler);
  Printf.printf "Launcher  : %s\n%!" (Launcher.to_string  launcher );


  (* Create the node file *)
  (*
  let () =
    let server_file =
      Filename.concat ezfio_filename "nodefile"
    in
    Out_channel.with_file server_file ~f:(fun out_channel ->
      Launcher.create_nodefile ()
      |> Out_channel.output_string out_channel
    )
  *)


  (* Get the configuration of executables *)
  let qmcchem =
    Lazy.force Qmcchem_config.qmcchem
  and qmc =
        [ Lazy.force Qmcchem_config.qmcchem ; "run" ; "-q" ]
  in


  if (start_dataserver) then
    begin
      (* Reset socket address in EZFIO *)
      Ezfio.set_simulation_http_server "tcp://127.0.0.1:65534";


      (* Start the data server *)
      let prog, args =
        qmcchem,  [| qmcchem; "run" ; "-d" ; ezfio_filename |]
      in
      let pid_dataserver =
        Watchdog.fork_exec ~prog ~args ()
      in
      Printf.printf "%7d : %s\n%!" pid_dataserver (String.concat " " (Array.to_list args))
    end;


  (* Check if the Zmq Rep socket is open *)
  let test_open_rep_socket () =
    let zmq_context =
      Zmq.Context.create ()
    in
    let socket =
      Zmq.Socket.create zmq_context Zmq.Socket.req
    and address =
      Ezfio.get_simulation_http_server ()
    in
    Zmq.Socket.set_receive_timeout socket 100;
    let reply =
       try
        (
           Zmq.Socket.connect socket address;
           Zmq.Socket.send socket (Message.(to_string Test));
           Zmq.Socket.recv socket
         ) with
       | Unix.Unix_error _ ->
           begin
             "Failed"
           end
    in
    Zmq.Socket.set_linger_period socket 1 ;
    Zmq.Socket.close socket;
    Zmq.Context.terminate zmq_context;
    reply = "OK"
  in


  (* Wait until the rep socket is open *)
  let rec count = function
  | 0  -> false
  | -1 -> true
  | n  ->
    if (not (test_open_rep_socket ())) then
      begin
       Unix.sleep 2;
       count (n-1);
      end
    else
      count (-1);
  in
  if (not (count 300)) then
    Watchdog.kill ();
  Unix.sleep 3;


  (* Start the qmc processes *)
  let prog, args_list =
    let launcher =
      Launcher.(find () |> to_string)
    in
    match launcher
    |> String.split_on_char ' '
    |> List.rev_map String.trim
    |> List.rev
    |> List.filter (fun x -> x <> "")
    with
    | launcher_exe::launcher_flags ->
       launcher_exe, launcher_exe :: launcher_flags @ qmc @ [
         Ezfio.get_simulation_http_server () ; ezfio_filename ]
    | _ -> failwith "Error in launcher"
  in
  let args = Array.of_list args_list in
  let pid_qmc =
    try
      Watchdog.fork_exec ~prog ~args ()
    with
    | Unix.Unix_error _ ->
        begin
          let command =
            String.concat " " args_list
          in
          Printf.printf "
============================================================
Error: Unable to run the following command
 %s
============================================================
\n%!" command ;
          Watchdog.kill ()
        end
  in
  Printf.printf "%7d : %s\n%!" pid_qmc (String.concat " " args_list);

  (* Wait for processes to finish *)
  Watchdog.join ()


let data_run ezfio_filename =
  Qmcchem_dataserver.run ezfio_filename ~daemon:false

let qmc_run dataserver ezfio_filename =
  Qmcchem_forwarder.run ezfio_filename dataserver

let ssh_run host dataserver ezfio_filename =
  print_endline ("ssh "^host^" "^ezfio_filename^" "^dataserver)

let run a d ?q ?s ezfio_filename =

  Qputils.set_ezfio_filename ezfio_filename;

  (* Signal handler to Kill properly all the processes *)
  let handler s =
    Printf.printf "QMC=Chem received signal %d... killing\n%!" s;
    Watchdog.kill ();
  in
  List.iter (fun s -> ignore @@ Sys.signal s (Sys.Signal_handle handler))
    [
      Sys.sigint  ;
      Sys.sigterm ;
      Sys.sigquit ;
    ]
  ;


  (* Validate input *)
  Input.validate ();
(*  Printf.printf "MD5 : %s\n" (Lazy.force Md5.hash) ; *)

  let runtype =
    match (a,d,q,s) with
    | (false,false, None, None) -> `Run
    | (false,true, None, None) -> `Data
    | (true,false, None, None) -> `Add
    | (false,false, Some dataserver, None) -> `Qmc dataserver
    | (false,false, Some dataserver, Some host) -> `Ssh (host, dataserver)
    | _ -> failwith "Options (-a|-d|-q [-s]) are mutually exclusive"
  in

  let run =
    match runtype with
    | `Run  -> full_run ~start_dataserver:true
    | `Data -> data_run
    | `Add  -> full_run ~start_dataserver:false
    | `Qmc dataserver -> qmc_run dataserver
    | `Ssh (host,dataserver)  -> ssh_run host dataserver
  in
  run ezfio_filename




let command () =
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - QMC=Chem command");
    set_description_doc "Run a calculation";

    [ { short='a' ; long="add" ; opt=Optional ;
        doc="Add more resources to a running calculation" ;
        arg=Without_arg ; };

      { short='d' ; long="data-server" ; opt=Optional ;
        doc="Start a dataserver process on the local host" ;
        arg=Without_arg ; };

      { short='q' ; long="local-qmc" ; opt=Optional ;
        doc="Start a qmc process on the local host attached to the addres given as an argument" ;
        arg=With_arg "<string>" ; };

      { short='s' ; long="remote-qmc" ; opt=Optional ;
        doc="Start a qmc process on the remote host as an argument" ;
        arg=With_arg "<string>" ; };

      anonymous "EZFIO_DIR" Mandatory "EZFIO directory";
    ]
    |> set_specs
  end;

  let a = Command_line.get_bool "add" in
  let d = Command_line.get_bool "data-server" in
  let q = Command_line.get "local-qmc" in
  let s = Command_line.get "remote-qmc" in

  let ezfio_file =
    match Command_line.anon_args () with
    | ezfio_file :: [] -> ezfio_file
    | _ -> (Command_line.help () ; failwith "Inconsistent command line")
  in

  run a d ?q ?s ezfio_file



