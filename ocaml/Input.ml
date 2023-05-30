open Qptypes
open Qputils



module Pseudo: sig

  type t = bool
  val doc   : string
  val read  : unit -> t
  val to_bool : t -> bool
  val of_bool : bool -> t
  val to_int  : t -> int
  val of_int  : int -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = bool

  let doc = "Compute pseudo-potentials"

  let of_bool x = x

  let to_bool x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_pseudo_do_pseudo ())) then
      Ezfio.set_pseudo_do_pseudo false;
    Ezfio.get_pseudo_do_pseudo ()
    |> of_bool


  let to_string t =
    to_bool t
    |> string_of_bool


  let of_string t =
    try
      String.lowercase_ascii t
      |> bool_of_string
      |> of_bool
    with
    | Invalid_argument msg -> failwith msg


  let to_int t =
    let t =
      to_bool t
    in
    if t then 1
    else 0


  let of_int = function
    | 0 -> false
    | 1 -> true
    | _ -> failwith "Expected 0 or 1"


end



module Fitcusp_factor : sig

  type t = float
  val doc   : string
  val read  : unit -> t
  val write : t -> unit
  val to_float : t -> float
  val of_float : float -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = float

  let doc = "Correct wave function to verify electron-nucleus cusp condition.
Fit is done for r < r_c(f) where r_c(f) = (1s orbital radius) x f. Value of f"

  let of_float x =
    if (x < 0.) then
       failwith "Fitcusp_factor should be >= 0.";
    if (x > 10.) then
       failwith "Fitcusp_factor is too large.";
    x

  let to_float x = x

  let read () =
    ignore @@
      Lazy.force Qputils.ezfio_filename ;
    if (not (Ezfio.has_simulation_nucl_fitcusp_factor ())) then
      begin
        let factor =
          Lazy.force Default.simulation_nucl_fitcusp_factor ;
        in
        Ezfio.set_simulation_nucl_fitcusp_factor factor
      end ;
    Ezfio.get_simulation_nucl_fitcusp_factor ()
    |> of_float


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_float t
    |> Ezfio.set_simulation_nucl_fitcusp_factor


  let to_string t =
    to_float t
    |> string_of_float


  let of_string t =
    try
      float_of_string t
      |> of_float
    with
    | Invalid_argument msg -> failwith msg

end

module Block_time : sig

  type t = int
  val doc   : string
  val read  : unit -> t
  val write : t -> unit
  val to_int : t -> int
  val of_int : int -> t
  val to_string : t -> string
  val of_string : string -> t
  val to_float  : t -> float
  val of_float  : float-> t

end = struct

  type t = int

  let doc = "Time (seconds) of a block"

  let of_int x =
    if (x < 1) then
      failwith "Block time should be >=1";
    if (x > 100000000) then
      failwith "Block time is too large (<= 100000000)";
    x


  let to_int x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_block_time ())) then
      Lazy.force Default.simulation_block_time
      |> Ezfio.set_simulation_block_time ;
    Ezfio.get_simulation_block_time ()
    |> of_int


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_int t
    |> Ezfio.set_simulation_block_time


  let to_string t =
    to_int t
    |> string_of_int


  let of_string t =
    int_of_string t
    |> of_int


  let to_float  t =
    to_int t
    |> float_of_int


  let of_float t =
    int_of_float t
    |> of_int


end

module Walk_num : sig

  type t = int
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_int : t -> int
  val of_int : int -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = int
  let doc = "Number of walkers per CPU core"

  let of_int x =
    if (x < 1) then
      failwith "Number of walkers should be >=1";
    if (x > 100_000) then
      failwith "Number of walkers is too large (<= 100_000)";
    x


  let to_int x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_electrons_elec_walk_num () )) then
      Lazy.force Default.electrons_elec_walk_num
      |> Ezfio.set_electrons_elec_walk_num ;
    Ezfio.get_electrons_elec_walk_num ()
    |> of_int


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_int t
    |> Ezfio.set_electrons_elec_walk_num


  let to_string t =
    to_int t
    |> string_of_int


  let of_string t =
    int_of_string t
    |> of_int


end

module Walk_num_tot : sig

  type t = int
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_int : t -> int
  val of_int : int -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = int
  let doc = "Total number of stored walkers for restart"

  let of_int x =
    if (x < 2) then
      failwith "Total number of stored walkers should be > 1";
    if (x > 100_000_000) then
      failwith "Number of walkers to store too large (<= 100.10^6)";
    x


  let to_int x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_electrons_elec_walk_num_tot () )) then
      Lazy.force Default.electrons_elec_walk_num_tot
      |> Ezfio.set_electrons_elec_walk_num_tot ;
    Ezfio.get_electrons_elec_walk_num_tot ()
    |> of_int


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_int t
    |> Ezfio.set_electrons_elec_walk_num_tot


  let to_string t =
    to_int t
    |> string_of_int


  let of_string t =
    int_of_string t
    |> of_int


end


module Stop_time : sig

  type t = int
  val read  : unit -> t
  val doc : string
  val write : t -> unit
  val to_int : t -> int
  val of_int : int -> t
  val to_float : t -> float
  val of_float : float -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = int

  let doc = "Requested simulation time (seconds)"

  let of_int x =
    if (x < 1) then
      failwith "Simulation time too short (>=1 s)";
    x


  let to_int x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_stop_time ())) then
      Lazy.force Default.simulation_stop_time
      |> Ezfio.set_simulation_stop_time ;
    Ezfio.get_simulation_stop_time ()
    |> of_int


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_int t
    |> Ezfio.set_simulation_stop_time


  let to_string t =
    to_int t
    |> string_of_int


  let of_string t =
    int_of_string t
    |> of_int


  let to_float  t =
    to_int t
    |> float_of_int


  let of_float t =
    int_of_float t
    |> of_int

end



module Method : sig

  type t = VMC | DMC | SRMC | FKMC | PDMC
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = VMC | DMC | SRMC | FKMC | PDMC

  let doc = "QMC Method : [ VMC | DMC | SRMC | FKMC | PDMC ]"

  let of_string = function
  | "VMC"  | "vmc"  -> VMC
  | "DMC"  | "dmc"  -> DMC
  | "SRMC" | "srmc" -> SRMC
  | "PDMC" | "pdmc" -> PDMC
  | "FKMC" | "fkmc" -> FKMC
  | x -> failwith ("Method should be [ VMC | DMC | SRMC | FKMC | PDMC ], not "^x^".")


  let to_string = function
  | VMC  -> "VMC"
  | DMC  -> "DMC"
  | SRMC -> "SRMC"
  | PDMC -> "PDMC"
  | FKMC -> "FKMC"


  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_method ())) then
      Lazy.force Default.simulation_method
      |> Ezfio.set_simulation_method ;
    Ezfio.get_simulation_method ()
    |> of_string


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_string t
    |> Ezfio.set_simulation_method


end



module Sampling : sig

  type t = Brownian | Langevin
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = Brownian | Langevin

  let doc = "Sampling algorithm : [ Langevin | Brownian ]"

  let of_string s =
    match String.capitalize_ascii (String.trim s) with
    | "Langevin" -> Langevin
    | "Brownian" -> Brownian
    | x -> failwith ("Sampling should be [ Brownian | Langevin ], not "^x^".")


  let to_string = function
  | Langevin -> "Langevin"
  | Brownian -> "Brownian"


  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_sampling ())) then
      Lazy.force Default.simulation_sampling
      |> Ezfio.set_simulation_sampling ;
    Ezfio.get_simulation_sampling ()
    |> of_string


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_string t
    |> Ezfio.set_simulation_sampling


end



module Trial_wf_energy : sig

  type t = float
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_float  : t -> float
  val of_float  : float -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = float
  let doc = "Energy of the trial wave function (au)"

  let of_float x =
    if (x > 0.) then
      failwith "Reference energy should not be positive.";
    if (x <= -1_000_000.) then
      failwith "Reference energy is too low.";
    x


  let to_float x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_e_trial ())) then
      to_float 0.
      |> Ezfio.set_simulation_e_trial;
    Ezfio.get_simulation_e_trial ()
    |> of_float


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_float t
    |> Ezfio.set_simulation_e_trial


  let of_string x =
     float_of_string x
     |> of_float


  let to_string x =
    to_float x
    |> string_of_float


end

module Ref_energy : sig

  type t = float
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_float  : t -> float
  val of_float  : float -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = float
  let doc = "Fixed reference energy to normalize DMC weights (au)"

  let of_float x =
    if (x > 0.) then
      failwith "Reference energy should not be positive.";
    if (x <= -1_000_000.) then
      failwith "Reference energy is too low.";
    x


  let to_float x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_e_ref ())) then
      to_float 0.
      |> Ezfio.set_simulation_e_ref;
    Ezfio.get_simulation_e_ref ()
    |> of_float


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_float t
    |> Ezfio.set_simulation_e_ref


  let of_string x =
     float_of_string x
     |> of_float


  let to_string x =
    to_float x
    |> string_of_float


end


module CI_threshold : sig

  type t = float
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_float  : t -> float
  val of_float  : float -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = float
  let doc = "Truncation t of the wave function : Remove determinants with a
contribution to the norm less than t (au)"

  let of_float x =
    if (x >= 1.) then
      failwith "Truncation of the wave function should be < 1.";
    if (x < 0.) then
      failwith "Truncation of the wave function should be positive.";
    x


  let to_float x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_ci_threshold ())) then
      Lazy.force Default.simulation_ci_threshold
      |> Ezfio.set_simulation_ci_threshold ;
    Ezfio.get_simulation_ci_threshold ()
    |> of_float


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_float t
    |> Ezfio.set_simulation_ci_threshold


  let of_string x =
    float_of_string x
    |> of_float


  let to_string x =
    to_float x
    |> string_of_float

end

module SRMC_projection_time : sig

  type t = float
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_float  : t -> float
  val of_float  : float -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = float
  let doc = "SRMC projection time (au)"

  let of_float x =
    if (x >= 100.) then
      failwith "SRMC Projection time should be < 100.";
    if (x <= 0.) then
      failwith "SRMC Projection time should be positive.";
    x


  let to_float x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_srmc_projection_time())) then
      Lazy.force Default.simulation_srmc_projection_time
      |> Ezfio.set_simulation_srmc_projection_time ;
    Ezfio.get_simulation_srmc_projection_time ()
    |> of_float


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_float t
    |> Ezfio.set_simulation_srmc_projection_time


  let of_string x =
    float_of_string x
    |> of_float


  let to_string x =
    to_float x
    |> string_of_float

end

module Time_step : sig

  type t = float
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_float  : t -> float
  val of_float  : float -> t
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = float
  let doc = "Simulation time step (au)"

  let of_float x =
    if (x >= 10.) then
      failwith "Time step should be < 10.";
    if (x <= 0.) then
      failwith "Time step should be positive.";
    x


  let to_float x = x

  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_simulation_time_step ())) then
      Lazy.force Default.simulation_time_step
      |> Ezfio.set_simulation_time_step ;
    Ezfio.get_simulation_time_step ()
    |> of_float


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_float t
    |> Ezfio.set_simulation_time_step


  let of_string x =
    float_of_string x
    |> of_float


  let to_string x =
    to_float x
    |> string_of_float

end

module Jastrow_type : sig

  type t = None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur
  let doc = "Type of Jastrow factor [ None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur ]"

  let of_string s =
    match String.capitalize_ascii (String.trim  s) with
    | "Core" -> Core
    | "Simple" -> Simple
    | "None" -> None
    | "Qmckl" -> Qmckl
    | "Mu" -> Mu
    | "Mu_1b" -> Mu_1b
    | "Muenv" -> Muenv
    | "Mur" -> Mur
    | _ -> failwith "Jastrow type should be [ None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur ]"


  let to_string = function
  | Core -> "Core"
  | Simple -> "Simple"
  | Mu -> "Mu"
  | Qmckl -> "Qmckl"
  | Mu_1b -> "Mu_1b"
  | Muenv -> "Muenv"
  | Mur -> "Mur"
  | None -> "None"


  let read () =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_jastrow_jast_type ())) then
      Lazy.force Default.jastrow_jast_type
      |> Ezfio.set_jastrow_jast_type ;
    Ezfio.get_jastrow_jast_type ();
    |> of_string


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_string t
    |> Ezfio.set_jastrow_jast_type


end

module Jpsi_type : sig

  type t = None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur
  let doc = "Type of Jpsi factor [ None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur ]"

  let of_string s = 
    match String.capitalize_ascii (String.trim  s) with
    | "Core" -> Core
    | "Simple" -> Simple
    | "None" -> None
    | "Qmckl" -> Qmckl
    | "Mu" -> Mu
    | "Mu_1b" -> Mu_1b
    | "Muenv" -> Muenv
    | "Mur" -> Mur
    | _ -> failwith "Jpsi type should be [ None | Core | Simple | Qmckl | Mu | Mu_1b | Muenv | Mur ]"


  let to_string = function
  | Core -> "Core"
  | Simple -> "Simple"
  | Mu -> "Mu"
  | Qmckl -> "Qmckl"
  | Mu_1b -> "Mu_1b"
  | Muenv -> "Muenv"
  | Mur -> "Mur"
  | None -> "None"


  let read () = 
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    if (not (Ezfio.has_jastrow_jpsi_type ())) then
      Lazy.force Default.jastrow_jpsi_type
      |> Ezfio.set_jastrow_jpsi_type ;
    Ezfio.get_jastrow_jpsi_type ();
    |> of_string


  let write t =
    let _ =
      Lazy.force Qputils.ezfio_filename
    in
    to_string t
    |> Ezfio.set_jastrow_jpsi_type


end

module Properties: sig

  type t = (Property.t * bool) list
  val doc : string
  val read  : unit -> t
  val write : t -> unit
  val to_string : t -> string
  val of_string : string -> t

end = struct

  type t = (Property.t * bool) list

  let doc =
    "Properties to sample. (X) is true and ( ) is false"


  let read () =
    List.rev_map (fun x -> (x, Property.calc x)) Property.all
    |> List.rev


  let write l =
    List.iter (fun (x,b) -> Property.set_calc x b) l


  let to_string l =
    List.rev_map (fun (x,b) ->
      let ch =
         if b then "X" else " "
      in
      Printf.sprintf "(%s) %s" ch (Property.to_string x)) l
    |> List.rev
    |> String.concat "\n"


  let of_string s =
    String.split_on_char '\n' s
    |> List.rev_map (fun x ->
       let (calc,prop) =
         String.trim x
         |> String_ext.rsplit2_exn ~on:' '
       in
       let prop =
         String.trim prop
         |> Property.of_string
       and calc =
         match calc with
         | "(X)" -> true
         | "( )" -> false
         | _ -> failwith " (X) or ( ) expected"
       in
       (prop, calc)
    )
    |> List.rev

end

(** Check if everything is correct in the input file. *)
let validate () =

  let _ =
    Lazy.force Qputils.ezfio_filename
  in

  (* Check if walkers are present *)
  if (not (Ezfio.has_electrons_elec_coord_pool ())) then
    Printf.printf "Warning: No initial walkers\n";

  let meth =
     Method.read ()
  and sampling =
     Sampling.read ()
  and ts =
     Time_step.read ()
  and do_pseudo =
    Pseudo.read ()
  in

  (* Check sampling and time steps *)
  let () =
    match (sampling, meth, Pseudo.to_bool do_pseudo) with
    | (Sampling.Brownian, Method.VMC, _) ->
      if ( (Time_step.to_float ts) >= 10. ) then
          warn "Time step seems large for VMC."
    | (Sampling.Langevin, Method.VMC, _) ->
      if ( (Time_step.to_float ts) <= 0.01 ) then
          warn "Time step seems small for Langevin sampling."
    | (Sampling.Brownian, _, true) ->
      if ( (Time_step.to_float ts) >= 0.5 ) then
          warn ( "Time step seems large for "^(Method.to_string meth) )
    | (Sampling.Brownian, _, false) ->
      if ( (Time_step.to_float ts) >= 0.01 ) then
          warn ( "Time step seems large for "^(Method.to_string meth) )
    | (Sampling.Langevin, _, _) ->
        failwith "Lanvegin sampling is incompatible with DMC"
  in


  (* Check E_ref is not zero *)
  let () =
    match (meth, Ref_energy.(read () |> to_float) ) with
    | (Method.SRMC,0.)
    | (Method.PDMC,0.)
    | (Method.FKMC,0.)
    | (Method.DMC,0.) -> failwith ("E_ref should not be zero in "^(Method.to_string meth) )
    | _          -> ()
  in

  (* Set block and total time*)
  let () =
    if ( (Block_time.read ()) > Stop_time.read ()) then
       failwith "Block time is longer than total time"
  in

  (* Check if E_loc if computed *)
  let () =
    match (meth, Property.(calc E_loc)) with
    | (Method.SRMC, false)
    | (Method.PDMC, false)
    | (Method.FKMC, false)
    | (Method.DMC, false) -> failwith ( "E_loc should be sampled in "^(Method.to_string meth) )
    | (Method.VMC, false) -> warn "Sampling of E_loc is not activated in input"
    | _ -> ()
  in

  (* Fitcusp is incompatible with pseudo *)
  let () =
    let f =
      Fitcusp_factor.read ()
      |> Fitcusp_factor.to_float
    in
    match (Pseudo.to_bool do_pseudo, f > 0.) with
    | (true, true) ->
      begin
        warn "Electron-nucleus cusp fitting is incompatible with Pseudopotentials.";
        Fitcusp_factor.of_float 0.
        |> Fitcusp_factor.write
      end
    | _ -> ()
  in


  (* Other Checks *)
  let () =
    let _ =
      Walk_num.read ()
    and _ =
      Walk_num_tot.read ()
    and _ =
      CI_threshold.read ()
    in ()
  in
  ()



