let split_re = 
  Str.regexp " +"


let split s =
  String.trim s
  |> Str.split split_re


let set_ezfio_filename ezfio_filename =
  let () = 
    if (not (Sys.file_exists ezfio_filename)) then
      failwith (ezfio_filename^" does not exist")
  in
  let () =
    if Sys.file_exists ezfio_filename && Sys.is_directory ezfio_filename then
      Ezfio.set_file ezfio_filename 
    else
      failwith ("Error : "^ezfio_filename^" is not a directory")
  in 
  let dir, result =
    Filename.dirname ezfio_filename,
    Filename.basename ezfio_filename
  in
  Unix.chdir dir ;
  Ezfio.set_file result
  

let ezfio_filename = lazy (
  let f =
    !Ezfio.ezfio_filename
  in
  let full_path = 
    match f with
    | "EZFIO_File" ->
        begin
          let args =
            Command_line.anon_args ()
            |> Array.of_list
          in
          if (Array.length args < 1) then
             failwith "Error : EZFIO directory not specified on the command line\n";
          args.(0)
        end
    | f -> f
  in
  set_ezfio_filename full_path;

  !Ezfio.ezfio_filename
)


let elec_num = lazy (
  Ezfio.set_file (Lazy.force ezfio_filename);
  Ezfio.get_electrons_elec_alpha_num () +
  Ezfio.get_electrons_elec_beta_num  ()
)


let walk_num = lazy (
  Ezfio.set_file (Lazy.force ezfio_filename);
  Ezfio.get_electrons_elec_walk_num ()
)


let warn msg =
  Printf.printf "Warning : %s\n%!" msg

let () =
  Random.self_init ()

