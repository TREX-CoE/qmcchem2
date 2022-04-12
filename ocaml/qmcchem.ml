let update_command_line () =
  let last = (Array.length Sys.argv) - 2 in
  Sys.argv.(0) <- Sys.argv.(0) ^ "_" ^ Sys.argv.(1);
  for i=1 to last do
    Sys.argv.(i) <- Sys.argv.(i+1)
  done;
  Sys.argv.(last+1) <- ""


let help () =
  Printf.printf "
qmcchem - QMC=Chem command

Usage:

  qmcchem [-h] COMMAND

Arguments:

  COMMAND       QMC=Chem command to run :
                [run|edit|stop|result|md5|info|debug]

Options:

  -h  --help    Prints the help message.

Description:

  Driver for subcommands.

"

let () =
  if Array.length Sys.argv < 2 then
    (help (); failwith "Inconsistent command line") ;

  match String.trim Sys.argv.(1) with
  | "-h" | "--help" ->
    begin
      help () ;
      exit 0
    end
  | _ ->
    begin
      let command =
        Sys.argv.(1)
      in
      update_command_line ();

      match command with
      | "debug"  -> let open Qmcchem_debug  in command ()
      | "edit"   -> let open Qmcchem_edit   in command ()
      | "info"   -> let open Qmcchem_info   in command ()
      | "md5"    -> let open Qmcchem_md5    in command ()
      | "result" -> let open Qmcchem_result in command ()
      | "run"    -> let open Qmcchem_run    in command ()
      | "stop"   -> let open Qmcchem_stop   in command ()
      | _        -> (help () ; failwith "Inconsistent command line")
    end

