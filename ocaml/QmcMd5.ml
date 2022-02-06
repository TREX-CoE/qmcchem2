(** Directory containing the list of input files. The directory is created is inexistant. *)
let input_directory = lazy (

  let ezfio_filename =
    Lazy.force Qputils.ezfio_filename
  in

  let dirname = 
    Filename.concat ezfio_filename "input"
  in

  begin
    if not (Sys.file_exists dirname) then 
      Unix.mkdir dirname 0o777
  end ;

  dirname
)


(** List of files responsible for the MD5 key of the input *)
let files_to_track = [
  "ao_basis/ao_coef.gz" ;
  "ao_basis/ao_expo.gz" ;
  "ao_basis/ao_nucl.gz" ;
  "ao_basis/ao_num" ; 
  "ao_basis/ao_power.gz" ;
  "ao_basis/ao_prim_num.gz" ; 
  "electrons/elec_alpha_num" ;
  "electrons/elec_beta_num" ;
  "electrons/elec_walk_num" ;
  "jastrow/jast_type" ;
  "mo_basis/mo_coef.gz" ;
  "mo_basis/mo_tot_num" ;
  "nuclei/nucl_charge.gz" ;
  "nuclei/nucl_coord.gz" ;
  "nuclei/nucl_num" ;
  "simulation/ci_threshold" ;
  "simulation/nucl_fitcusp_factor" ;
  "simulation/jast_a_up_dn" ;
  "simulation/jast_a_up_up" ;
  "simulation/jast_b_up_dn" ;
  "simulation/jast_b_up_up" ;
  "simulation/jast_core_a1" ;
  "simulation/jast_core_a2" ;
  "simulation/jast_core_b1" ;
  "simulation/jast_core_b2" ;
  "simulation/jast_een_e_a.gz" ;
  "simulation/jast_een_e_b.gz" ;
  "simulation/jast_een_n.gz" ;
  "simulation/jast_pen.gz" ;
  "simulation/method" ;
  "simulation/time_step" ;
  "simulation/dmc_projection_time" ;
  "spindeterminants/bit_kind" ;
  "spindeterminants/n_det" ;
  "spindeterminants/n_det_alpha" ;
  "spindeterminants/n_det_beta" ;
  "spindeterminants/n_int" ;
  "spindeterminants/n_states" ;
  "spindeterminants/psi_coef_matrix_columns.gz" ;
  "spindeterminants/psi_coef_matrix_rows.gz" ;
  "spindeterminants/psi_coef_matrix_values.gz" ;
  "spindeterminants/psi_det_alpha.gz" ;
  "spindeterminants/psi_det_beta.gz" ;
  "/pseudo/do_pseudo" ;
  "/pseudo/mo_pseudo_grid.gz" ;
  "/pseudo/pseudo_dz_kl.gz";
  "/pseudo/pseudo_klocmax" ;
  "/pseudo/pseudo_n_k.gz" ;
  "/pseudo/pseudo_v_kl.gz" ;
  "/pseudo/pseudo_grid_rmax" ;
  "/pseudo/pseudo_kmax" ;
  "/pseudo/pseudo_n_kl.gz" ;
  "/pseudo/pseudo_dz_k.gz" ;
  "/pseudo/pseudo_grid_size" ;
  "/pseudo/pseudo_v_k.gz" ;
  ]


(** Get an MD5 ke from the content of a file. *)
let hash_file filename =
  if Sys.file_exists filename then
    begin
      let ic = open_in filename in
      let result = 
        Cryptokit.hash_channel (Cryptokit.Hash.md5 ()) ic
        |> Cryptokit.transform_string (Cryptokit.Hexa.encode ()) 
      in
      close_in ic;
      result
    end
  else ""
  

(** Cache containing the current value of the MD5 hash. *)
let _hash = 
  ref None

(** Get the hash correcponding to the EZFIO file. *)
let hash () = 
  let compute_hash () =
    let ezfio_filename =
      Lazy.force Qputils.ezfio_filename
    in
    let old_md5 =
      if Ezfio.has_simulation_md5_key () then
         Ezfio.get_simulation_md5_key ()
      else
         ""
    in
    let md5_string = 
      files_to_track
      |> List.rev_map (fun x -> Printf.sprintf "%s/%s" ezfio_filename x)
      |> List.rev_map hash_file
      |> String.concat ""
    in

    let new_md5 = 
      md5_string
      |> Cryptokit.hash_string (Cryptokit.Hash.md5 ()) 
      |> Cryptokit.transform_string (Cryptokit.Hexa.encode ()) 
    in
    if (new_md5 <> old_md5) then
      begin
        Printf.eprintf "Info : MD5 key changed\n   %s\n-> %s\n%!" old_md5 new_md5 ;
        Ezfio.set_simulation_md5_key new_md5
      end
    ;
    new_md5
  in
  match (!_hash) with
  | Some key -> key
  | None -> 
      begin
        let key =
          compute_hash ()
        in
        _hash := Some key ;
        key
      end

(** Reset the cache of the MD5 hash. *)
let reset_hash () =
  _hash := None;
  ignore (hash ())



