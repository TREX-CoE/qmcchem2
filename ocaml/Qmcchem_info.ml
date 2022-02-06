

let run ezfio_filename = 
  Qputils.set_ezfio_filename ezfio_filename;
  let qmcchem_info =
      Lazy.force Qmcchem_config.qmcchem_info
  in
  let prog, argv = 
      qmcchem_info, 
    [| qmcchem_info ; ezfio_filename |]
  in
  ignore @@
    Unix.execvp prog argv 

let command () =
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - QMC=Chem command");
    set_description_doc "Display info on an EZFIO database";
    [ anonymous "EZFIO_DIR" Mandatory "EZFIO directory" ]
    |> set_specs
  end;

  let ezfio_file =
    match Command_line.anon_args () with
    | ezfio_file :: [] -> ezfio_file                                                          
    | _ -> (Command_line.help () ; failwith "Inconsistent command line")
  in

  run ezfio_file

