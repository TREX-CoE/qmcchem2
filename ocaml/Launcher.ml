type t =
  | Srun
  | MPI
  | Bash

let to_string = function
  | Srun     -> String.concat " " [ "srun" ; 
                    try Sys.getenv "QMCCHEM_SRUN_FLAGS"  with Not_found -> ""
                  ]
  | Bash     -> "env"
  | MPI      -> String.concat " " [ Lazy.force Qmcchem_config.mpirun ;
                    try Sys.getenv "QMCCHEM_MPIRUN_FLAGS"  with Not_found -> ""
                  ]

(*
let to_string = function
  | MPI
  | Srun     -> String.concat ~sep:" " [ Lazy.force Qmcchem_config.mpirun ;
      match Sys.getenv "QMCCHEM_MPIRUN_FLAGS" with
      | None -> ""
      | Some p -> p
    ]
  | Bash     -> "env"
  *)



(** Find the launcher for the current job scheduler *)
let find () =

  let result =
    match Scheduler.find () with
    | Scheduler.SLURM -> Srun
    | Scheduler.Batch
    | Scheduler.PBS
    | Scheduler.SGE ->
      if Lazy.force Qmcchem_config.has_mpirun then
        MPI
      else
        Bash
  in
  result


(** Create a file contaning the list of nodes and the number of available CPUs *)
let create_nodefile () =

  let launcher =
    find ()
  in

  let launcher_command =
    to_string launcher
  in

  let h =
    Hashtbl.create 1000
  in

  let in_channel =
    Unix.open_process_in (launcher_command^" hostname -s")
  in
  String_ext.input_lines in_channel
  |> List.map  String.trim
  |> List.iter ( fun host ->
       let n =
          match Hashtbl.find_opt h host with
          | Some x -> x+1
          | None -> 1
       in
       Hashtbl.replace h host n
     );
  match
    Unix.close_process_in in_channel
  with
  | _ -> ();


  let f =
    match launcher with
    | MPI ->
        fun (node, n) ->
            Printf.sprintf "%s slots=%d\n" node n
    | Srun
    | Bash ->
        fun (node, n) ->
            Printf.sprintf "%s %d\n" node n
  in
  Hashtbl.fold (fun k v a -> (f (k,v)) :: a) h []
  |> String.concat "\n"




