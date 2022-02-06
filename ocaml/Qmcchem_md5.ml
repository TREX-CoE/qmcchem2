
let run ?c ?d ~l ~update ezfio_filename =

  Qputils.set_ezfio_filename ezfio_filename;

  let input_directory =
    Lazy.force QmcMd5.input_directory
  in

  let handle_options () =

    let current_md5 = 
      QmcMd5.hash ()
    in

    let filename_of_key key = 
      Filename.concat  input_directory  key
    in

    let key_is_valid key =
      let filename = 
        filename_of_key key
      in
      Sys.file_exists filename
    in

    if (update) then
       begin
          Printf.printf "Updating\n%!" ;
          let update_one old_key =
            Qmcchem_edit.run ~c:false ~input:(filename_of_key old_key) ezfio_filename;
            QmcMd5.reset_hash ();
            let new_key =
              QmcMd5.hash ()
            in

            if (old_key <> new_key) then
              begin
                let prefix = 
                  Filename.concat ezfio_filename "blocks" 
                in
                let new_name = 
                  Filename.concat prefix new_key 
                and old_name = 
                  Filename.concat prefix old_key 
                in
                Printf.printf "Renaming %s -> %s\n" old_name new_name;
                try Sys.rename old_name new_name with
                | Sys_error _ -> ();

                let old_name = 
                  String.concat "/" [ ezfio_filename; "input"; old_key ]
                in
                Printf.printf "Removing %s\n%!" old_name;
                try Sys.remove old_name with
                | Sys_error _ -> ();
            end
          in
          Sys.readdir input_directory
          |> Array.iter (fun x -> update_one x) ;
          Printf.printf "Done\n%!" ;
       end
    ;

    let () =
      match c with
      | None -> ()
      | Some new_md5 -> 
          if (key_is_valid new_md5) then
            Qmcchem_edit.run ~c:false ~input:(filename_of_key new_md5) ezfio_filename
          else
            failwith ("Error: " ^ new_md5 ^ " does not exist") 
    in

    let () = 
      match l with
      | false -> ()
      | true  ->
          Sys.readdir   input_directory
          |> Array.iter (fun md5 -> 
              let filename =
                Filename.concat  input_directory  md5
              in
              let this =
                if (md5 = current_md5) then
                  "<-"
                else
                  ""
              in
              let date =
                let open Unix in
                localtime (stat filename).st_mtime 
                |> Time.string_of_date
              in
              Printf.printf "%s : %s  %s\n" md5 date this)
    in

    let () =
      match d with
      | None -> ()
      | Some other_key ->
          if (key_is_valid other_key) then
            let command =
              String.concat " "
               [ "diff" ; "-u" ; "-w" ;
                 (filename_of_key current_md5) ;
                 (filename_of_key other_key) ]
            in
            ignore @@ Unix.system command
          else
            failwith ("Error: " ^ other_key ^ " does not exist") 
    in
    ()

  in

  match (c,d,l,update) with
  | (None,None,false,false) -> 
      Printf.printf "Current key :\n%s\n" (QmcMd5.hash ())
  | _ -> handle_options ()


let command () =
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - QMC=Chem command");
    set_description_doc "Manipulate input MD5 keys";

    [ { short='c' ; long="clear" ; opt=Optional ;
	doc="Change to input to <key>" ;
        arg=With_arg "<string>" ; };

      { short='d' ; long="diff" ; opt=Optional ;
	doc="Show input differences with <key>" ;
        arg=With_arg "<string>" ; };

      { short='l' ; long="list" ; opt=Optional ;
	doc="List all the saved MD5 keys." ;
        arg=Without_arg ; };

      { short='u' ; long="update" ; opt=Optional ;
	doc="Update to the latest MD5 format." ;
        arg=Without_arg ; };

      anonymous "EZFIO_DIR" Mandatory "EZFIO directory";
    ]
    |> set_specs 
  end;

  let update = Command_line.get_bool "update" in
  let c = Command_line.get "clear" in
  let d = Command_line.get "diff" in
  let l = Command_line.get_bool "list" in

  let ezfio_file =
    match Command_line.anon_args () with
    | ezfio_file :: [] -> ezfio_file
    | _ -> (Command_line.help () ; failwith "Inconsistent command line")
  in

  run ?c ?d ~l ~update ezfio_file 




