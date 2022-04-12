let run ~t ezfio_filename=

  Qputils.set_ezfio_filename ezfio_filename;

  if (not (Ezfio.has_simulation_http_server ())) then
     failwith "QMC=Chem is not running"
  ;

  let zmq_context =
    Zmq.Context.create ()
  in

  Printf.printf "Debugging %s\n%!" ezfio_filename;
  let socket =
    Zmq.Socket.create zmq_context Zmq.Socket.sub
  in

  let address =
    match (Ezfio.get_simulation_http_server ()
        |> String_ext.rsplit2 ~on:':' )
    with
    | Some (a,p) -> a^":"^( (int_of_string p)+4 |> string_of_int )
    | None -> failwith "Badly formed address"
  in
  Zmq.Socket.connect socket address;
  Zmq.Socket.subscribe socket "";

  if t then
      begin
        let re_split =
           Str.regexp " *: *"
        in
        let tot_size =
          ref 0.
        in
        while true
        do
          let msg =
            Zmq.Socket.recv socket
          in
          let (socket, bytes)  =
            match Str.split re_split msg with
            |  socket :: bytes :: _ ->
                (socket, float_of_string bytes)
            | _ -> (print_endline msg ; ("", 0.))
          in
          tot_size := !tot_size +. bytes;
          Printf.printf "%f\n%!" !tot_size;
          Unix.sleep 1
        done
      end
  else
      begin
        while true
        do
          let msg =
            Zmq.Socket.recv socket
          in
          Printf.printf "%s\n%!" msg;
        done
      end



let command () =
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - QMC=Chem command");
    set_description_doc "Debug ZeroMQ communications";
    [ { short='t' ; long="traffic" ; opt=Optional ;
        doc="Print traffic in bytes" ;
        arg=Without_arg } ;

        anonymous "EZFIO_DIR" Mandatory "EZFIO directory" ]
    |> set_specs
  end;

  let t = Command_line.get_bool "traffic" in
  let ezfio_file =
    match Command_line.anon_args () with
    | ezfio_file :: [] -> ezfio_file
    | _ -> (Command_line.help () ; failwith "Inconsistent command line")
  in

  run t ezfio_file





