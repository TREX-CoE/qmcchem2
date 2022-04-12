let bind_socket ~socket_type ~socket ~address =
  let rec loop = function
    | 0 -> failwith @@ Printf.sprintf
        "Unable to bind the forwarder's %s socket : %s\n"
        socket_type address
    | -1 -> ()
    | i ->
      try
        Zmq.Socket.bind socket address;
        loop (-1)
      with
      | Unix.Unix_error _ -> (Unix.sleep 1 ; loop (i-1) )
      | other_exception -> raise other_exception
  in loop 10



let run ezfio_filename dataserver =

  let dataserver_address, dataserver_port =
    String.sub dataserver 6 (String.length dataserver - 6)
    |> String_ext.lsplit2_exn ~on:':'
  and qmc =
    Lazy.force Qmcchem_config.qmc
  in

  (* Go into /dev/shm *)
  Unix.chdir Qmcchem_config.dev_shm;

  let tmpdir =
    ezfio_filename ^ "_" ^ dataserver_port
  in

  (* Port of the data server *)
  let port =
    (int_of_string dataserver_port)+10
  in

  (* Build qmc executable command *)
  let prog, argv =
      qmc,
    [| qmc ; ezfio_filename ;
       Printf.sprintf "ipc://%s:%d" Qmcchem_config.dev_shm port |];
  in

  (* Create the temporary directory. If it is possible, then the process is a
   * master and the forwarder will start. Otherwise, only start a qmc process.
   *)
  let () =
    try
      Unix.mkdir tmpdir 0o755;
      Unix.chdir tmpdir
    with
    | Unix.Unix_error _ ->
      begin
        Unix.chdir tmpdir;
        Unix.sleep 1;
        if Sys.file_exists "PID" then
          begin
            let pid =
              let ic = open_in "PID" in
              try
                int_of_string (input_line ic)
              with
              | End_of_file -> -1
            in
            match pid with
            | -1 -> ()
            | pid ->
                try
                  Unix.kill pid 0 ;
                  ignore @@ Unix.execvp prog argv
                with
                | Unix.Unix_error (Unix.ESRCH, _, _) -> ()
          end
      end
  in

  (* Now, only one forwarder will execute the following code *)
  let oc = open_out "PID" in
  Unix.getpid ()
  |> Printf.sprintf "%d\n"
  |> output_string oc
  ; close_out oc;

  (* Fork a qmc *)
  ignore @@
    Watchdog.fork_exec ~prog ~args:argv ();

  (* Fetch input *)
  let zmq_context =
    Zmq.Context.create ()
  in

  let terminate () =
    (* Clean up the temp directory *)
    Unix.chdir Qmcchem_config.dev_shm;
    Zmq.Context.terminate zmq_context ;
    for i=port to port+4
    do
      let filename =
         Printf.sprintf ":%d" i
      in
      try
         Sys.remove filename
      with
      | _ ->  ()
      ;
    done;
    let command =
      Printf.sprintf "rm -rf -- \"%s\" " tmpdir
    in
    try
      ignore @@ Unix.system command
    with
    | Unix.Unix_error _  -> print_endline "Unable to remove temporary directory"
    ;
    Watchdog.kill ()
  in


  (* Signal handler to Kill properly all the processes *)
  let handler s =
    Printf.printf "Forwarder received signal %d... killing\n%!" s;
    terminate ()
  in
  List.iter (fun s -> ignore @@ Sys.signal s (Sys.Signal_handle handler))
    [
      Sys.sigint  ;
      Sys.sigterm ;
      Sys.sigquit ;
    ]
  ;


  (* Fetch walkers *)
  let walk_num =
    ref 0
  and walkers =
    ref []
  in


  (* Status thread *)
  let status =
    ref Status.Running
  in

  let start_status_thread =

    let f () =
      let pub_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.pub
      and address =
        Printf.sprintf "ipc://%s:%d" Qmcchem_config.dev_shm (port+1);
      in
      bind_socket "PUB" pub_socket address;

      let sub_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.sub
      and address =
        Printf.sprintf "tcp://%s:%d" dataserver_address (port+1-10)
      in
      Zmq.Socket.connect sub_socket address;
      Zmq.Socket.subscribe sub_socket "";

      let pollitem =
        Zmq.Poll.mask_of
        [| (sub_socket, Zmq.Poll.In) ;
        |]
      in

      while (!status <> Status.Stopped)
      do
        let polling =
          Zmq.Poll.poll ~timeout:1000 pollitem
        in
        if (polling.(0) = Some Zmq.Poll.In) then
          begin
            let msg =
              Zmq.Socket.recv ~block:false sub_socket
            in
            Zmq.Socket.send pub_socket msg;
            status := Status.of_string msg;
          end;
      done;
      List.iter (fun socket ->
        Zmq.Socket.set_linger_period socket 1000 ;
        Zmq.Socket.close socket)
        [ sub_socket ; pub_socket ]
    in
    Thread.create f
  in

  let start_log_thread =

    let f () =
      let sub_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.xsub
      and address =
        Printf.sprintf "ipc://%s:%d" Qmcchem_config.dev_shm (port+3);
      in
      bind_socket "XSUB" sub_socket address;

      let pub_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.xpub
      and address =
        Printf.sprintf "tcp://%s:%d" dataserver_address (port+3-10)
      in
      Zmq.Socket.connect pub_socket address;

      let pollitem =
        Zmq.Poll.mask_of
        [| (sub_socket, Zmq.Poll.In) ;
           (pub_socket, Zmq.Poll.In) ;
        |]
      in

      (* Main loop *)
      while (!status <> Status.Stopped)
      do
        let polling =
          Zmq.Poll.poll ~timeout:1000 pollitem
        in
        if (polling.(0) = Some Zmq.Poll.In) then
          begin
            Zmq.Socket.recv ~block:false sub_socket
            |> Zmq.Socket.send pub_socket ;
          end
        else if (polling.(1) = Some Zmq.Poll.In) then
          begin
            Printf.eprintf "Forwarder subscribe\n%!";
            Zmq.Socket.recv ~block:false pub_socket
            |> Zmq.Socket.send sub_socket ;
          end
      done;
      List.iter (fun socket ->
        Zmq.Socket.set_linger_period socket 1000 ;
        Zmq.Socket.close socket)
        [ sub_socket ; pub_socket ]
    in
    Thread.create f
  in

  (* Proxy thread *)
  let start_proxy_thread =
    let f () =

      let req_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.req
      in
      Zmq.Socket.connect req_socket dataserver;
      Zmq.Socket.set_receive_timeout req_socket 600_000;

      let dealer_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.dealer
      in

      bind_socket "PROXY" dealer_socket "inproc://dealer";
      Zmq.Socket.set_receive_high_water_mark dealer_socket 100_000;
      Zmq.Socket.set_send_high_water_mark dealer_socket 100_000;
      Zmq.Socket.set_immediate dealer_socket true;
      Zmq.Socket.set_linger_period dealer_socket 600_000;

      let fetch_walkers () =
        Zmq.Socket.send_all req_socket ["get_walkers" ; string_of_int !walk_num ];
        Zmq.Socket.recv_all req_socket
      in

      let pollitem =
        Zmq.Poll.mask_of
        [| (dealer_socket, Zmq.Poll.In) ;
        |]
      in

      (* EZFIO Cache *)
      let ezfio_cache =
        Hashtbl.create 63
      in
      let handle_ezfio msg =
        match Hashtbl.find_opt ezfio_cache msg with
        | Some result -> result
        | None ->
          begin
            Zmq.Socket.send_all req_socket ["Ezfio" ; msg];
            let result =
              Zmq.Socket.recv_all req_socket
            in
            Hashtbl.add ezfio_cache msg result;
            result
          end
      in


      (* Main loop *)
      while (!status <> Status.Stopped)
      do
        let polling =
          Zmq.Poll.poll ~timeout:1000 pollitem
        in
        if (polling.(0) = Some Zmq.Poll.In) then
          begin
            let raw_msg =
              Zmq.Socket.recv_all ~block:false dealer_socket
            in
            let header, msg =
              let rec aux header = function
                | ""   :: msg  -> List.rev ("" :: header), Message.create msg
                | head :: tail -> aux (head::header) tail
                | _            -> failwith "Too many routers in the middle"
              in
              aux [] (List.rev @@ List.rev_map String.trim raw_msg)
            in
            let handle message =
              match message with
              | Message.Ezfio ezfio_msg ->
                  let result =
                    handle_ezfio ezfio_msg
                  in
                  Zmq.Socket.send_all dealer_socket (header @ result)
              | Message.GetWalkers n_walks ->
                begin
                  if (!walk_num = 0) then
                    begin
                      walk_num := Qptypes.Strictly_positive_int.to_int n_walks;
                      walkers := fetch_walkers ();
                    end;
                  Zmq.Socket.send_all dealer_socket (header @ !walkers);
                  walkers := fetch_walkers ();
                end
              | Message.Test ->
                  Zmq.Socket.send_all dealer_socket (header @ [ "OK" ])
              | Message.Error _ ->  ()
              | Message.Register _
              | Message.Unregister _
              | Message.Walkers  _
              | Message.Property _ ->
                  failwith "Bad message"
            in handle msg
          end;
      done;
      Zmq.Socket.set_linger_period dealer_socket 1000 ;
      Zmq.Socket.set_linger_period req_socket 1000 ;
      Zmq.Socket.close dealer_socket;
      Zmq.Socket.close req_socket;
    in
    Thread.create f
  in

  (* Main thread *)
  let start_main_thread =
    let f () =

      let dealer_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.dealer
      in
      Zmq.Socket.connect dealer_socket dataserver;
      Zmq.Socket.set_linger_period dealer_socket 600_000;

      let proxy_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.dealer
      in
      Zmq.Socket.connect proxy_socket "inproc://dealer";

      let router_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.router
      and address =
        Printf.sprintf "ipc://%s:%d" Qmcchem_config.dev_shm (port);
      in
      bind_socket "ROUTER" router_socket address;
      Zmq.Socket.set_receive_high_water_mark router_socket 100000;
      Zmq.Socket.set_send_high_water_mark router_socket 100000;
      Zmq.Socket.set_immediate router_socket true;
      Zmq.Socket.set_linger_period router_socket 600_000;

      (* Pull socket for computed data *)
      let push_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.push
      and address =
        Printf.sprintf "tcp://%s:%d" dataserver_address (port+2-10)
      in
      Zmq.Socket.connect push_socket address;
      Zmq.Socket.set_linger_period push_socket 600_000;

      let pull_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.pull
      and address =
        Printf.sprintf "ipc://%s:%d" Qmcchem_config.dev_shm (port+2);
      in
      bind_socket "PULL" pull_socket address;


      (* Handles messages coming into the ROUTER socket. *)
      let handle_router () =
        let raw_msg =
          Zmq.Socket.recv_all ~block:false router_socket
        in
        let header, msg =
          let rec aux header = function
            | ""   :: msg  -> List.rev ("" :: header), Message.create msg
            | head :: tail -> aux (head::header) tail
            | _            -> failwith "Too many routers in the middle"
          in
          aux [] (List.rev @@ List.rev_map String.trim raw_msg)
        in
        let handle message =
          match message with
          | Message.GetWalkers _
          | Message.Ezfio _
          | Message.Test ->
              Zmq.Socket.send_all proxy_socket raw_msg
          | Message.Register _
          | Message.Unregister _ ->
              Zmq.Socket.send_all dealer_socket raw_msg
          | Message.Walkers (_, _, _)
          | Message.Property _ ->
              failwith "Bad message"
          | Message.Error _ -> ()
        in handle msg
      in

      let handle_dealer () =
        Zmq.Socket.recv_all ~block:false dealer_socket
        |> Zmq.Socket.send_all router_socket
      in

      let handle_proxy () =
        Zmq.Socket.recv_all ~block:false proxy_socket
        |> Zmq.Socket.send_all router_socket
      in

      let select_n_of ~n ~len l =
        let a =
          Array.of_list l
        in
        let s =
          (Array.length a)/ len
        in
        let fetch i =
          let rec loop accu = function
          | -1 -> accu
          | k -> loop ((Array.get a (i*len+k)) :: accu) (k-1)
          in
          loop [] (len-1)
        in
        let rec select accu = function
        | 0 -> accu
        | i -> let new_accu =
                (fetch @@ Random.int s) :: accu
              in
              select new_accu (i-1)
        in
        select [] n
        |> List.concat
      in

      (* Handles messages coming into the PULL socket. *)
      let handle_pull () =
        let message =
          Zmq.Socket.recv_all ~block:false pull_socket
        in
        let new_message =
          match message with
          | "elec_coord":: hostname :: pid :: id :: n_str :: rest ->
            let n =
              int_of_string n_str
            in
            let len =
              if !walk_num = 0 then n else
              n / !walk_num
            in
            if (n < 5*len) then
              message
            else
              List.concat [ [ "elec_coord" ; hostname ; pid ; id ;
              string_of_int (5*len)] ; ( select_n_of ~n:5 ~len rest ) ]
          | prop :: c :: pid :: b :: d :: w :: [] -> message
          | prop :: c :: pid :: b :: d :: w :: l ->
            if Qmcchem_config.binary_io then
              match Message.create message with
              | Message.Property block ->
                  prop :: c :: pid :: b :: d :: w :: "bin" ::
                  (Block.to_bytes block |> Bytes.unsafe_to_string ) :: []
              | _ -> failwith "Inconsistent message"
            else
              message
          | _ -> message
        in
        Zmq.Socket.send_all  push_socket  new_message
      in

      (* Polling item to poll ROUTER and PULL sockets. *)
      let pollitem =
        Zmq.Poll.mask_of
        [| (router_socket , Zmq.Poll.In) ;
            (pull_socket  , Zmq.Poll.In) ;
            (dealer_socket, Zmq.Poll.In) ;
            (proxy_socket , Zmq.Poll.In)
        |]
      in
      (* Main loop *)
      while (!status <> Status.Stopped)
      do
        let polling =
          Zmq.Poll.poll ~timeout:1000 pollitem
        in
        if (polling.(0) = Some Zmq.Poll.In) then
          handle_router ();
        if (polling.(1) = Some Zmq.Poll.In) then
          handle_pull ();
        if (polling.(2) = Some Zmq.Poll.In) then
          handle_dealer ();
        if (polling.(3) = Some Zmq.Poll.In) then
          handle_proxy ();
      done;
      List.iter (fun socket ->
        Zmq.Socket.set_linger_period socket 1000 ;
        Zmq.Socket.close socket)
      [ router_socket ; dealer_socket ; push_socket ; pull_socket ; proxy_socket ]
    in
    Thread.create f
  in


  (* Start the status thread and the main thread *)
  begin
    try
      (List.iter Thread.join
        [ start_status_thread ();
          start_log_thread ();
          start_proxy_thread ();
          start_main_thread ();
        ])
    with
    | err ->
      begin
        print_endline "Trapped error. Waiting 10 seconds...";
        status := Status.Stopping;
        Unix.sleep 10 ;
        raise err
      end
  end;

  (* Wait for the qmc process to complete *)
  try
    ignore (Watchdog.join ());
    terminate ()
  with
  | error ->
    begin
      terminate ();
      raise error
    end

