open Qputils

type t =
| Stopped
| Queued
| Running
| Stopping
;;

let of_string = function
| "Stopped"   -> Stopped
| "Queued"    -> Queued
| "Running"   -> Running
| "Stopping"  -> Stopping
| _           -> failwith "Invalid status"
;;

let of_int    = function
| 0 -> Stopped
| 1 -> Queued
| 2 -> Running
| 3 -> Stopping
| _ -> failwith "Invalid status"
;;

let to_string = function
| Stopped   -> "Stopped"
| Queued    -> "Queued"
| Running   -> "Running"
| Stopping  -> "Stopping"
;;

let to_int = function
| Stopped   -> 0
| Queued    -> 1
| Running   -> 2
| Stopping  -> 3
;;

let is_set = ref false

let read () =
  if not !is_set then
    (Ezfio.set_file (Lazy.force ezfio_filename); is_set := true);
  Ezfio.get_simulation_do_run ()
  |> of_int
;;

let write x =
  if not !is_set then
    (Ezfio.set_file (Lazy.force ezfio_filename); is_set := true);
  to_int x
  |> Ezfio.set_simulation_do_run
;;

