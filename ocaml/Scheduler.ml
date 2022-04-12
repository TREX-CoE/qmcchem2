type t =
  | SGE
  | PBS
  | SLURM
  | Batch


let to_string = function
  | SGE   -> "SGE"
  | PBS   -> "PBS"
  | SLURM -> "SLURM"
  | Batch -> "Batch"



let find () =
  let scheduler =
    [  "SLURM_NODELIST" ; "PE_HOSTFILE" ; "PBS_NODEFILE" ]
    |> List.map (function x ->
         try ignore @@ (Sys.getenv x) ; Some x with
         | Not_found -> None
         )
    |> List.hd
  in
  let result =
    match scheduler with
      | Some "SLURM_NODELIST"  -> SLURM
      | Some "PE_HOSTFILE"     -> SGE
      | Some "PBS_NODEFILE"    -> PBS
      | None                   -> Batch
      | Some x  -> failwith (Printf.sprintf "Scheduler %s not found" x)
  in
  result



