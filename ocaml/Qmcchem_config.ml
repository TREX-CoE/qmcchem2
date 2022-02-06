
(** QMC=Chem installation directory *)
let root = lazy (
  try Sys.getenv "QMCCHEM_PATH" with
  | Not_found -> failwith "QMCCHEM_PATH environment variable not set"
)


(* PATH environment variable as a list of strings *)
let path =  lazy (
  let p =
    try Sys.getenv "PATH" with
    | Not_found -> failwith "PATH environment variable is not set"
  in
  String.split_on_char ':' p
)



(* Full path of a binary taken from the PATH *)
let full_path exe =
  let rec in_path_rec = function
  | [] -> None
  | head :: tail ->
    begin
      let fp =
        Filename.concat head exe
      in
      if Sys.file_exists fp then
        Some fp
      else
        in_path_rec tail
    end
  in
  Lazy.force path
  |> in_path_rec



(* True if an executable is in the PATH *)
let in_path x =
   match full_path x with
   | Some _ -> true
   | None   -> false


let has_parallel   = lazy( in_path "parallel" )
let has_mpirun     = lazy( in_path "mpirun"   )
let has_srun       = lazy( in_path "parallel" )
let has_qmc        = lazy( in_path "qmc"      )


let mpirun = lazy (
    try
      Sys.getenv "QMCCHEM_MPIRUN"
    with
    | Not_found -> "mpirun"
)

let qmcchem = lazy(
  Filename.concat (Lazy.force root) "bin/qmcchem"
)
and qmc = lazy(
  Filename.concat (Lazy.force root) "bin/qmc"
)
and qmcchem_info = lazy(
  Filename.concat (Lazy.force root) "bin/qmcchem_info"
)

and qmc_create_walkers = lazy(
  Filename.concat (Lazy.force root) "bin/qmc_create_walkers"
)

let dev_shm = "/dev/shm/"

(** Name of the host on which the data server runs *)
let hostname = lazy (
  try
    Unix.gethostname ()
  with
  | _ -> "127.0.0.1"
)


external get_ipv4_address_for_interface : string -> string =
  "get_ipv4_address_for_interface"


let ip_address = lazy (
  let interface =
    try Some (Sys.getenv "QMCCHEM_NIC")
    with Not_found -> None
  in
  match interface with
  | None ->
      begin
        try
          let host =
            Lazy.force hostname
            |> Unix.gethostbyname
          in
          Unix.string_of_inet_addr host.h_addr_list.(0);
        with
        | Unix.Unix_error _ ->
            failwith "Unable to find IP address from host name."
      end
  | Some interface ->
        let result = get_ipv4_address_for_interface interface in
        if String.sub result 0 5 = "error" then
          Printf.sprintf "Unable to use network interface %s" interface
          |> failwith
        else
          result
)



let binary_io =
  try
    let qmcchem_io = Sys.getenv "QMCCHEM_IO" in
    qmcchem_io = "B" || qmcchem_io = "b"
  with Not_found -> false

