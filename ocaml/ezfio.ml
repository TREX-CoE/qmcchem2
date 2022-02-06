(*
    EZFIO is an automatic generator of I/O libraries
    Copyright (C) 2009 Anthony SCEMAMA, CNRS

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    Anthony Scemama
    LCPQ - IRSAMC - CNRS
    Universite Paul Sabatier
    118, route de Narbonne
    31062 Toulouse Cedex 4
    scemama@irsamc.ups-tlse.fr
*)

let version = "2.0.7";;
let library = "/home/scemama/qmcchem/EZFIO";;

(*
Exceptions
==========
*)

exception Read_only of string

(*
State variables
===============
*)
let read_only = ref false
let ezfio_filename = ref "EZFIO_File"

(*
Helper functions
=================
*)

let check_readonly =
  if !read_only then
    raise (Read_only "Read only attribute is set")


let exists path =
  let filename = Filename.concat path ".version" in
  Sys.file_exists filename


let has group name =
  let dirname = Filename.concat !ezfio_filename group in
  if (exists dirname) then
    Sys.file_exists @@ Filename.concat dirname name
  else
    false


let has_array group name = has group (name^".gz")

let mkdir path =
  check_readonly;
  if (exists path) then
    raise (Failure (path^" already exists"));
  Unix.mkdir path 0o777;
  let out_channel = open_out (path^"/.version") in
  Printf.fprintf out_channel "%s\n" version ;
  close_out out_channel


let create_group group =
  let dirname = Filename.concat !ezfio_filename group in
  if not (exists dirname) then
      mkdir dirname



let maxval l =
  match l with
  | [] -> None
  | [a] -> Some a
  | hd::tail -> Some (List.fold_left max hd tail)


let minval l =
  match l with
  | [] -> None
  | [a] -> Some a
  | hd::tail -> Some (List.fold_left min hd tail)


let at arr idx = arr.(idx)

(*
let size (_) = 0
*)

let n_count_ch (str,_,v) =
  let rec do_work accu = function
  | 0 -> accu
  | i ->
    let newaccu =
      if str.[i-1] == v then accu+1
      else accu
    in do_work newaccu (i-1)
  in do_work 0 (String.length str)


let n_count (l,_,v) =
  let rec do_work accu = function
  | []  -> accu
  | h::tail ->
    let newaccu =
      if h == v then accu+1
      else accu
    in do_work newaccu tail
  in do_work 0 l


(*
Scalars
=======
*)

(*
Read
----
*)

let read_scalar type_conversion group name =
  let in_filename = Filename.concat !ezfio_filename @@ Filename.concat group name in
  let in_channel = open_in in_filename in
  let trimmed_line = String.trim (input_line in_channel) in
  let result = type_conversion trimmed_line in
  begin
    close_in in_channel ;
    result
  end


let fortran_bool_of_string = function
  | "T" | "t" -> true
  | "F" | "f" -> false
  | x -> raise (Failure ("fortran_bool_of_string should be T or F: "^x))


let fortran_string_of_bool = function
 | true -> "T\n"
 | false-> "F\n"


let read_int   = read_scalar   int_of_string
let read_int64 = read_scalar Int64.of_string
let read_float = read_scalar float_of_string
let read_string= read_scalar (fun (x:string) -> x)
let read_bool  = read_scalar fortran_bool_of_string

(*
Write
-----
*)

let print_int    out_channel i = Printf.fprintf out_channel "%20d\n" i
let print_int64  out_channel i = Printf.fprintf out_channel "%20Ld\n" i
let print_float  out_channel f = Printf.fprintf out_channel "%24.15e\n" f
let print_string out_channel s = Printf.fprintf out_channel "%s\n" s
let print_bool   out_channel b = Printf.fprintf out_channel "%s\n" (fortran_string_of_bool b)

let write_scalar print_fun group name s =
  check_readonly;
  create_group group;
  let out_filename = Filename.concat !ezfio_filename @@ Filename.concat group name in
  let out_channel = open_out out_filename in
  begin
     print_fun out_channel s;
     close_out out_channel
  end


let write_int    = write_scalar print_int
let write_int64  = write_scalar print_int64
let write_float  = write_scalar print_float
let write_bool   = write_scalar print_bool
let write_string = write_scalar print_string




(*
Arrays
======
*)

type 'a ezfio_data =
| Ezfio_item of 'a array
| Ezfio_data of ('a ezfio_data) array



type 'a ezfio_array =
{ rank : int ;
  dim  : int array;
  data : 'a ezfio_data ;
}

let ezfio_array_of_array ~rank ~dim ~data =
  assert (rank > 0);
  let read_1d data nmax =
    (Ezfio_item (Array.sub data 0 nmax), Array.sub data nmax ((Array.length data) - nmax))
  in
  let rec read_nd data = function
    | m when m<1 -> raise (Failure "dimension should not be <1")
    | 1 -> read_1d data dim.(0)
    | m ->
      let rec do_work accu data = function
      | 0 -> (Array.of_list (List.rev accu), data)
      | n ->
        let (newlist,rest) = read_nd data (m-1) in
        do_work (newlist::accu) rest (n-1)
      in
      let (data,rest) = do_work [] data dim.(m-1) in
      (Ezfio_data data,rest)
  in
  let (result,_) = read_nd data rank in
  { rank= rank;
    dim=  dim;
    data=result;
  }


let ezfio_array_of_list ~rank ~dim ~data =
  assert (rank > 0);
  let read_1d data nmax =
    let rec do_work accu data = function
    | 0 -> (Array.of_list (List.rev accu), data)
    | n ->
      begin
        match data with
        | x::rest -> do_work (x::accu) rest (n-1)
        | [] -> raise (Failure "Array is not consistent")
      end
    in
    let (data,rest) = do_work [] data nmax in
    (Ezfio_item data,rest)
  in
  let rec read_nd data = function
    | m when m<1 -> raise (Failure "dimension should not be <1")
    | 1 -> read_1d data dim.(0)
    | m ->
      let rec do_work accu data = function
      | 0 -> (Array.of_list (List.rev accu), data)
      | n ->
        let (newlist,rest) = read_nd data (m-1) in
        do_work (newlist::accu) rest (n-1)
      in
      let (data,rest) = do_work [] data dim.(m-1) in
      (Ezfio_data data,rest)
  in
  let (result,_) = read_nd data rank in
  { rank= rank;
    dim=  dim;
    data=result;
  }



let ezfio_get_element { rank=r ; dim=d ; data=data } coord =
  (*assert ((List.length coord) == r);*)
  let rec do_work buffer = function
  | [c] ->
    begin match buffer with
    | Ezfio_item buffer -> buffer.(c)
    | Ezfio_data buffer -> raise (Failure "Error in ezfio_get_element")
    end
  | c::tail ->
    begin match buffer with
    | Ezfio_item buffer -> raise (Failure "Error in ezfio_get_element")
    | Ezfio_data buffer -> do_work buffer.(c) tail
    end
  | [] -> raise (Failure "Error in ezfio_get_element")
  in
  do_work data coord



let flattened_ezfio { rank ; dim ; data } =
  let flatten_2 d =
    let l = List.rev_map (function
      | Ezfio_item i -> i
      | Ezfio_data i -> assert false
    ) (List.rev (Array.to_list d))
    in Array.concat l
  in

  let rec flatten_n rank d =
    if (rank = 2) then
      flatten_2 d
    else
      let l = List.rev_map (function
        | Ezfio_data x -> flatten_n (rank-1) x
        | Ezfio_item _ -> assert false
      ) (List.rev (Array.to_list d))
      in Array.concat l
  in

  match data with
    | Ezfio_item d -> d
    | Ezfio_data d -> flatten_n rank d



(*
Read
----
*)

let read_rank in_channel =
    let trimmed_line = String.trim (input_line in_channel) in
    int_of_string trimmed_line


let read_dimensions in_channel =
    let line = input_line in_channel in
    let arr_of_str  =
      String.split_on_char ' ' line
      |> List.filter (fun x -> x <> "")
      |> Array.of_list
    in
    Array.map int_of_string arr_of_str



let read_array type_conversion group name : 'a ezfio_array =
  let in_filename = (Filename.concat !ezfio_filename @@ Filename.concat group name) ^ ".gz" in
  let in_channel = Unix.open_process_in ("zcat "^in_filename) in
  (* Read rank *)
  let rank = read_rank in_channel
  (* Read dimensions *)
  and dimensions = read_dimensions in_channel
  in
  begin
      assert (rank == Array.length dimensions) ;
      (* Read one-dimensional arrays *)
      let read_1d nmax =
        Ezfio_item (Array.init nmax (fun i->
          let buffer = String.trim (input_line in_channel) in
          try type_conversion buffer with
          | Failure s -> failwith (s^": "^buffer)))
      in
      (* Read multi-dimensional arrays *)
      let rec read_nd = function
        | m when m<1 -> raise (Failure "dimension should not be <1")
        | 1 -> read_1d dimensions.(0)
        | m ->
          let rec do_work accu = function
          | 0 -> Array.of_list (List.rev accu)
          | n ->
            let newlist = read_nd (m-1) in
            do_work (newlist::accu) (n-1)
          in
          Ezfio_data (do_work [] dimensions.(m-1))
      in
      let result = {
         rank = rank ;
         dim = dimensions ;
         data = read_nd rank ;
      }
      in
      match Unix.close_process_in in_channel with
        | Unix.WEXITED _ -> result
        | _ -> failwith ("Failed in reading compressed file "^in_filename)
  end


let read_int_array    = read_array          int_of_string
let read_int64_array  = read_array        Int64.of_string
let read_float_array  = read_array        float_of_string
let read_bool_array   = read_array fortran_bool_of_string
let read_string_array = read_array (fun (x:string) -> x)


(*
Write
-----
*)

let write_array print_fun group name a =
  check_readonly;
  create_group group;
  let out_filename = (Filename.concat !ezfio_filename @@ Filename.concat group name)^".gz" in
  let command = "gzip > "^out_filename in
  let out_channel = Unix.open_process_out command in
  let { rank=rank ; dim=dimensions ; data=data } = a in
  let data = flattened_ezfio a in
  begin
     (* Write rank *)
     Printf.fprintf out_channel "%3d\n" rank;
     (* Write dimensions *)
     Array.iter (Printf.fprintf out_channel "%20d ") dimensions;
     Printf.fprintf out_channel "\n%!";
     Array.iter (print_fun out_channel) data;
     match Unix.close_process_out out_channel with
      | Unix.WEXITED _ -> ()
      | _ -> failwith ("Failed writing compressed file "^out_filename)
  end


let write_int_array    = write_array print_int
let write_int64_array  = write_array print_int64
let write_float_array  = write_array print_float
let write_string_array = write_array print_string
let write_bool_array   = write_array print_bool

(*
Library routines
*)

let set_file filename =
  if not (exists filename) then
  begin
    mkdir filename;
    mkdir (Filename.concat filename "ezfio");
    let command = Printf.sprintf "
      LANG= date > %s/ezfio/creation
      echo $USER > %s/ezfio/user
      echo %s > %s/ezfio/library" filename filename library filename
    in
    if (Sys.command command <> 0) then
      raise (Failure ("Unable to create new ezfio file:\n"^filename))
  end ;
  ezfio_filename := filename



let get_ao_basis_ao_num () = read_int "ao_basis" "ao_num";;
let set_ao_basis_ao_num = write_int "ao_basis" "ao_num" ;;
let has_ao_basis_ao_num () = has "ao_basis" "ao_num" ;;
let get_ao_basis_ao_prim_num () = read_int_array "ao_basis" "ao_prim_num";;
let set_ao_basis_ao_prim_num = write_int_array "ao_basis" "ao_prim_num";;
let has_ao_basis_ao_prim_num () = has_array "ao_basis" "ao_prim_num" ;;
let get_ao_basis_ao_nucl () = read_int_array "ao_basis" "ao_nucl";;
let set_ao_basis_ao_nucl = write_int_array "ao_basis" "ao_nucl";;
let has_ao_basis_ao_nucl () = has_array "ao_basis" "ao_nucl" ;;
let get_ao_basis_ao_power () = read_int_array "ao_basis" "ao_power";;
let set_ao_basis_ao_power = write_int_array "ao_basis" "ao_power";;
let has_ao_basis_ao_power () = has_array "ao_basis" "ao_power" ;;
let get_ao_basis_ao_coef () = read_float_array "ao_basis" "ao_coef";;
let set_ao_basis_ao_coef = write_float_array "ao_basis" "ao_coef";;
let has_ao_basis_ao_coef () = has_array "ao_basis" "ao_coef" ;;
let get_ao_basis_ao_expo () = read_float_array "ao_basis" "ao_expo";;
let set_ao_basis_ao_expo = write_float_array "ao_basis" "ao_expo";;
let has_ao_basis_ao_expo () = has_array "ao_basis" "ao_expo" ;;
let get_mo_basis_mo_num () = read_int "mo_basis" "mo_num";;
let set_mo_basis_mo_num = write_int "mo_basis" "mo_num" ;;
let has_mo_basis_mo_num () = has "mo_basis" "mo_num" ;;
let get_mo_basis_mo_coef () = read_float_array "mo_basis" "mo_coef";;
let set_mo_basis_mo_coef = write_float_array "mo_basis" "mo_coef";;
let has_mo_basis_mo_coef () = has_array "mo_basis" "mo_coef" ;;
let get_mo_basis_mo_classif () = read_string_array "mo_basis" "mo_classif";;
let set_mo_basis_mo_classif = write_string_array "mo_basis" "mo_classif";;
let has_mo_basis_mo_classif () = has_array "mo_basis" "mo_classif" ;;
let get_mo_basis_mo_energy () = read_float_array "mo_basis" "mo_energy";;
let set_mo_basis_mo_energy = write_float_array "mo_basis" "mo_energy";;
let has_mo_basis_mo_energy () = has_array "mo_basis" "mo_energy" ;;
let get_mo_basis_mo_occ () = read_float_array "mo_basis" "mo_occ";;
let set_mo_basis_mo_occ = write_float_array "mo_basis" "mo_occ";;
let has_mo_basis_mo_occ () = has_array "mo_basis" "mo_occ" ;;
let get_mo_basis_mo_symmetry () = read_string_array "mo_basis" "mo_symmetry";;
let set_mo_basis_mo_symmetry = write_string_array "mo_basis" "mo_symmetry";;
let has_mo_basis_mo_symmetry () = has_array "mo_basis" "mo_symmetry" ;;
let get_electrons_elec_alpha_num () = read_int "electrons" "elec_alpha_num";;
let set_electrons_elec_alpha_num = write_int "electrons" "elec_alpha_num" ;;
let has_electrons_elec_alpha_num () = has "electrons" "elec_alpha_num" ;;
let get_electrons_elec_beta_num () = read_int "electrons" "elec_beta_num";;
let set_electrons_elec_beta_num = write_int "electrons" "elec_beta_num" ;;
let has_electrons_elec_beta_num () = has "electrons" "elec_beta_num" ;;
let get_electrons_elec_walk_num_tot () = read_int "electrons" "elec_walk_num_tot";;
let set_electrons_elec_walk_num_tot = write_int "electrons" "elec_walk_num_tot" ;;
let has_electrons_elec_walk_num_tot () = has "electrons" "elec_walk_num_tot" ;;
let get_electrons_elec_walk_num () = read_int "electrons" "elec_walk_num";;
let set_electrons_elec_walk_num = write_int "electrons" "elec_walk_num" ;;
let has_electrons_elec_walk_num () = has "electrons" "elec_walk_num" ;;
let get_electrons_elec_coord_pool () = read_float_array "electrons" "elec_coord_pool";;
let set_electrons_elec_coord_pool = write_float_array "electrons" "elec_coord_pool";;
let has_electrons_elec_coord_pool () = has_array "electrons" "elec_coord_pool" ;;
let get_electrons_elec_coord_pool_size () = read_int "electrons" "elec_coord_pool_size";;
let set_electrons_elec_coord_pool_size = write_int "electrons" "elec_coord_pool_size" ;;
let has_electrons_elec_coord_pool_size () = has "electrons" "elec_coord_pool_size" ;;
let get_electrons_elec_fitcusp_radius () = read_float "electrons" "elec_fitcusp_radius";;
let set_electrons_elec_fitcusp_radius = write_float "electrons" "elec_fitcusp_radius" ;;
let has_electrons_elec_fitcusp_radius () = has "electrons" "elec_fitcusp_radius" ;;
let get_nuclei_nucl_num () = read_int "nuclei" "nucl_num";;
let set_nuclei_nucl_num = write_int "nuclei" "nucl_num" ;;
let has_nuclei_nucl_num () = has "nuclei" "nucl_num" ;;
let get_nuclei_nucl_label () = read_string_array "nuclei" "nucl_label";;
let set_nuclei_nucl_label = write_string_array "nuclei" "nucl_label";;
let has_nuclei_nucl_label () = has_array "nuclei" "nucl_label" ;;
let get_nuclei_nucl_charge () = read_float_array "nuclei" "nucl_charge";;
let set_nuclei_nucl_charge = write_float_array "nuclei" "nucl_charge";;
let has_nuclei_nucl_charge () = has_array "nuclei" "nucl_charge" ;;
let get_nuclei_nucl_coord () = read_float_array "nuclei" "nucl_coord";;
let set_nuclei_nucl_coord = write_float_array "nuclei" "nucl_coord";;
let has_nuclei_nucl_coord () = has_array "nuclei" "nucl_coord" ;;
let get_spindeterminants_n_det_alpha () = read_int "spindeterminants" "n_det_alpha";;
let set_spindeterminants_n_det_alpha = write_int "spindeterminants" "n_det_alpha" ;;
let has_spindeterminants_n_det_alpha () = has "spindeterminants" "n_det_alpha" ;;
let get_spindeterminants_n_det_beta () = read_int "spindeterminants" "n_det_beta";;
let set_spindeterminants_n_det_beta = write_int "spindeterminants" "n_det_beta" ;;
let has_spindeterminants_n_det_beta () = has "spindeterminants" "n_det_beta" ;;
let get_spindeterminants_n_det () = read_int "spindeterminants" "n_det";;
let set_spindeterminants_n_det = write_int "spindeterminants" "n_det" ;;
let has_spindeterminants_n_det () = has "spindeterminants" "n_det" ;;
let get_spindeterminants_n_int () = read_int "spindeterminants" "n_int";;
let set_spindeterminants_n_int = write_int "spindeterminants" "n_int" ;;
let has_spindeterminants_n_int () = has "spindeterminants" "n_int" ;;
let get_spindeterminants_bit_kind () = read_int "spindeterminants" "bit_kind";;
let set_spindeterminants_bit_kind = write_int "spindeterminants" "bit_kind" ;;
let has_spindeterminants_bit_kind () = has "spindeterminants" "bit_kind" ;;
let get_spindeterminants_n_states () = read_int "spindeterminants" "n_states";;
let set_spindeterminants_n_states = write_int "spindeterminants" "n_states" ;;
let has_spindeterminants_n_states () = has "spindeterminants" "n_states" ;;
let get_spindeterminants_psi_det_alpha () = read_int64_array "spindeterminants" "psi_det_alpha";;
let set_spindeterminants_psi_det_alpha = write_int64_array "spindeterminants" "psi_det_alpha";;
let has_spindeterminants_psi_det_alpha () = has_array "spindeterminants" "psi_det_alpha" ;;
let get_spindeterminants_psi_det_beta () = read_int64_array "spindeterminants" "psi_det_beta";;
let set_spindeterminants_psi_det_beta = write_int64_array "spindeterminants" "psi_det_beta";;
let has_spindeterminants_psi_det_beta () = has_array "spindeterminants" "psi_det_beta" ;;
let get_spindeterminants_psi_coef_matrix_rows () = read_int_array "spindeterminants" "psi_coef_matrix_rows";;
let set_spindeterminants_psi_coef_matrix_rows = write_int_array "spindeterminants" "psi_coef_matrix_rows";;
let has_spindeterminants_psi_coef_matrix_rows () = has_array "spindeterminants" "psi_coef_matrix_rows" ;;
let get_spindeterminants_psi_coef_matrix_columns () = read_int_array "spindeterminants" "psi_coef_matrix_columns";;
let set_spindeterminants_psi_coef_matrix_columns = write_int_array "spindeterminants" "psi_coef_matrix_columns";;
let has_spindeterminants_psi_coef_matrix_columns () = has_array "spindeterminants" "psi_coef_matrix_columns" ;;
let get_spindeterminants_psi_coef_matrix_values () = read_float_array "spindeterminants" "psi_coef_matrix_values";;
let set_spindeterminants_psi_coef_matrix_values = write_float_array "spindeterminants" "psi_coef_matrix_values";;
let has_spindeterminants_psi_coef_matrix_values () = has_array "spindeterminants" "psi_coef_matrix_values" ;;
let get_spindeterminants_n_svd_coefs_unique () = read_int "spindeterminants" "n_svd_coefs_unique";;
let set_spindeterminants_n_svd_coefs_unique = write_int "spindeterminants" "n_svd_coefs_unique" ;;
let has_spindeterminants_n_svd_coefs_unique () = has "spindeterminants" "n_svd_coefs_unique" ;;
let get_spindeterminants_n_svd_coefs () = read_int "spindeterminants" "n_svd_coefs";;
let set_spindeterminants_n_svd_coefs = write_int "spindeterminants" "n_svd_coefs" ;;
let has_spindeterminants_n_svd_coefs () = has "spindeterminants" "n_svd_coefs" ;;
let get_spindeterminants_n_svd_selected () = read_int "spindeterminants" "n_svd_selected";;
let set_spindeterminants_n_svd_selected = write_int "spindeterminants" "n_svd_selected" ;;
let has_spindeterminants_n_svd_selected () = has "spindeterminants" "n_svd_selected" ;;
let get_spindeterminants_n_svd_toselect () = read_int "spindeterminants" "n_svd_toselect";;
let set_spindeterminants_n_svd_toselect = write_int "spindeterminants" "n_svd_toselect" ;;
let has_spindeterminants_n_svd_toselect () = has "spindeterminants" "n_svd_toselect" ;;
let get_spindeterminants_psi_svd_alpha_unique () = read_float_array "spindeterminants" "psi_svd_alpha_unique";;
let set_spindeterminants_psi_svd_alpha_unique = write_float_array "spindeterminants" "psi_svd_alpha_unique";;
let has_spindeterminants_psi_svd_alpha_unique () = has_array "spindeterminants" "psi_svd_alpha_unique" ;;
let get_spindeterminants_psi_svd_beta_unique () = read_float_array "spindeterminants" "psi_svd_beta_unique";;
let set_spindeterminants_psi_svd_beta_unique = write_float_array "spindeterminants" "psi_svd_beta_unique";;
let has_spindeterminants_psi_svd_beta_unique () = has_array "spindeterminants" "psi_svd_beta_unique" ;;
let get_spindeterminants_psi_svd_coefs_unique () = read_float_array "spindeterminants" "psi_svd_coefs_unique";;
let set_spindeterminants_psi_svd_coefs_unique = write_float_array "spindeterminants" "psi_svd_coefs_unique";;
let has_spindeterminants_psi_svd_coefs_unique () = has_array "spindeterminants" "psi_svd_coefs_unique" ;;
let get_spindeterminants_psi_svd_alpha () = read_float_array "spindeterminants" "psi_svd_alpha";;
let set_spindeterminants_psi_svd_alpha = write_float_array "spindeterminants" "psi_svd_alpha";;
let has_spindeterminants_psi_svd_alpha () = has_array "spindeterminants" "psi_svd_alpha" ;;
let get_spindeterminants_psi_svd_beta () = read_float_array "spindeterminants" "psi_svd_beta";;
let set_spindeterminants_psi_svd_beta = write_float_array "spindeterminants" "psi_svd_beta";;
let has_spindeterminants_psi_svd_beta () = has_array "spindeterminants" "psi_svd_beta" ;;
let get_spindeterminants_psi_svd_coefs () = read_float_array "spindeterminants" "psi_svd_coefs";;
let set_spindeterminants_psi_svd_coefs = write_float_array "spindeterminants" "psi_svd_coefs";;
let has_spindeterminants_psi_svd_coefs () = has_array "spindeterminants" "psi_svd_coefs" ;;
let get_spindeterminants_psi_svd_alpha_numselected () = read_int_array "spindeterminants" "psi_svd_alpha_numselected";;
let set_spindeterminants_psi_svd_alpha_numselected = write_int_array "spindeterminants" "psi_svd_alpha_numselected";;
let has_spindeterminants_psi_svd_alpha_numselected () = has_array "spindeterminants" "psi_svd_alpha_numselected" ;;
let get_spindeterminants_psi_svd_beta_numselected () = read_int_array "spindeterminants" "psi_svd_beta_numselected";;
let set_spindeterminants_psi_svd_beta_numselected = write_int_array "spindeterminants" "psi_svd_beta_numselected";;
let has_spindeterminants_psi_svd_beta_numselected () = has_array "spindeterminants" "psi_svd_beta_numselected" ;;
let get_spindeterminants_psi_svd_alpha_numtoselect () = read_int_array "spindeterminants" "psi_svd_alpha_numtoselect";;
let set_spindeterminants_psi_svd_alpha_numtoselect = write_int_array "spindeterminants" "psi_svd_alpha_numtoselect";;
let has_spindeterminants_psi_svd_alpha_numtoselect () = has_array "spindeterminants" "psi_svd_alpha_numtoselect" ;;
let get_spindeterminants_psi_svd_beta_numtoselect () = read_int_array "spindeterminants" "psi_svd_beta_numtoselect";;
let set_spindeterminants_psi_svd_beta_numtoselect = write_int_array "spindeterminants" "psi_svd_beta_numtoselect";;
let has_spindeterminants_psi_svd_beta_numtoselect () = has_array "spindeterminants" "psi_svd_beta_numtoselect" ;;
let get_simulation_do_run () = read_int "simulation" "do_run";;
let set_simulation_do_run = write_int "simulation" "do_run" ;;
let has_simulation_do_run () = has "simulation" "do_run" ;;
let get_simulation_stop_time () = read_int "simulation" "stop_time";;
let set_simulation_stop_time = write_int "simulation" "stop_time" ;;
let has_simulation_stop_time () = has "simulation" "stop_time" ;;
let get_simulation_equilibration () = read_bool "simulation" "equilibration";;
let set_simulation_equilibration = write_bool "simulation" "equilibration" ;;
let has_simulation_equilibration () = has "simulation" "equilibration" ;;
let get_simulation_http_server () = read_string "simulation" "http_server";;
let set_simulation_http_server = write_string "simulation" "http_server" ;;
let has_simulation_http_server () = has "simulation" "http_server" ;;
let get_simulation_do_jast () = read_bool "simulation" "do_jast";;
let set_simulation_do_jast = write_bool "simulation" "do_jast" ;;
let has_simulation_do_jast () = has "simulation" "do_jast" ;;
let get_simulation_nucl_fitcusp_factor () = read_float "simulation" "nucl_fitcusp_factor";;
let set_simulation_nucl_fitcusp_factor = write_float "simulation" "nucl_fitcusp_factor" ;;
let has_simulation_nucl_fitcusp_factor () = has "simulation" "nucl_fitcusp_factor" ;;
let get_simulation_method () = read_string "simulation" "method";;
let set_simulation_method = write_string "simulation" "method" ;;
let has_simulation_method () = has "simulation" "method" ;;
let get_simulation_block_time () = read_int "simulation" "block_time";;
let set_simulation_block_time = write_int "simulation" "block_time" ;;
let has_simulation_block_time () = has "simulation" "block_time" ;;
let get_simulation_sampling () = read_string "simulation" "sampling";;
let set_simulation_sampling = write_string "simulation" "sampling" ;;
let has_simulation_sampling () = has "simulation" "sampling" ;;
let get_simulation_save_data () = read_bool "simulation" "save_data";;
let set_simulation_save_data = write_bool "simulation" "save_data" ;;
let has_simulation_save_data () = has "simulation" "save_data" ;;
let get_simulation_time_step () = read_float "simulation" "time_step";;
let set_simulation_time_step = write_float "simulation" "time_step" ;;
let has_simulation_time_step () = has "simulation" "time_step" ;;
let get_simulation_print_level () = read_int "simulation" "print_level";;
let set_simulation_print_level = write_int "simulation" "print_level" ;;
let has_simulation_print_level () = has "simulation" "print_level" ;;
let get_simulation_ci_threshold () = read_float "simulation" "ci_threshold";;
let set_simulation_ci_threshold = write_float "simulation" "ci_threshold" ;;
let has_simulation_ci_threshold () = has "simulation" "ci_threshold" ;;
let get_simulation_md5_key () = read_string "simulation" "md5_key";;
let set_simulation_md5_key = write_string "simulation" "md5_key" ;;
let has_simulation_md5_key () = has "simulation" "md5_key" ;;
let get_simulation_e_ref () = read_float "simulation" "e_ref";;
let set_simulation_e_ref = write_float "simulation" "e_ref" ;;
let has_simulation_e_ref () = has "simulation" "e_ref" ;;
let get_simulation_e_trial () = read_float "simulation" "e_trial";;
let set_simulation_e_trial = write_float "simulation" "e_trial" ;;
let has_simulation_e_trial () = has "simulation" "e_trial" ;;
let get_simulation_srmc_projection_time () = read_float "simulation" "srmc_projection_time";;
let set_simulation_srmc_projection_time = write_float "simulation" "srmc_projection_time" ;;
let has_simulation_srmc_projection_time () = has "simulation" "srmc_projection_time" ;;
let get_jastrow_jast_type () = read_string "jastrow" "jast_type";;
let set_jastrow_jast_type = write_string "jastrow" "jast_type" ;;
let has_jastrow_jast_type () = has "jastrow" "jast_type" ;;
let get_jastrow_jast_a_up_up () = read_float "jastrow" "jast_a_up_up";;
let set_jastrow_jast_a_up_up = write_float "jastrow" "jast_a_up_up" ;;
let has_jastrow_jast_a_up_up () = has "jastrow" "jast_a_up_up" ;;
let get_jastrow_jast_a_up_dn () = read_float "jastrow" "jast_a_up_dn";;
let set_jastrow_jast_a_up_dn = write_float "jastrow" "jast_a_up_dn" ;;
let has_jastrow_jast_a_up_dn () = has "jastrow" "jast_a_up_dn" ;;
let get_jastrow_jast_b_up_up () = read_float "jastrow" "jast_b_up_up";;
let set_jastrow_jast_b_up_up = write_float "jastrow" "jast_b_up_up" ;;
let has_jastrow_jast_b_up_up () = has "jastrow" "jast_b_up_up" ;;
let get_jastrow_jast_b_up_dn () = read_float "jastrow" "jast_b_up_dn";;
let set_jastrow_jast_b_up_dn = write_float "jastrow" "jast_b_up_dn" ;;
let has_jastrow_jast_b_up_dn () = has "jastrow" "jast_b_up_dn" ;;
let get_jastrow_mu_erf () = read_float "jastrow" "mu_erf";;
let set_jastrow_mu_erf = write_float "jastrow" "mu_erf" ;;
let has_jastrow_mu_erf () = has "jastrow" "mu_erf" ;;
let get_jastrow_jast_pen () = read_float_array "jastrow" "jast_pen";;
let set_jastrow_jast_pen = write_float_array "jastrow" "jast_pen";;
let has_jastrow_jast_pen () = has_array "jastrow" "jast_pen" ;;
let get_jastrow_jast_een_e_a () = read_float_array "jastrow" "jast_een_e_a";;
let set_jastrow_jast_een_e_a = write_float_array "jastrow" "jast_een_e_a";;
let has_jastrow_jast_een_e_a () = has_array "jastrow" "jast_een_e_a" ;;
let get_jastrow_jast_een_e_b () = read_float_array "jastrow" "jast_een_e_b";;
let set_jastrow_jast_een_e_b = write_float_array "jastrow" "jast_een_e_b";;
let has_jastrow_jast_een_e_b () = has_array "jastrow" "jast_een_e_b" ;;
let get_jastrow_jast_een_n () = read_float_array "jastrow" "jast_een_n";;
let set_jastrow_jast_een_n = write_float_array "jastrow" "jast_een_n";;
let has_jastrow_jast_een_n () = has_array "jastrow" "jast_een_n" ;;
let get_jastrow_jast_core_a1 () = read_float_array "jastrow" "jast_core_a1";;
let set_jastrow_jast_core_a1 = write_float_array "jastrow" "jast_core_a1";;
let has_jastrow_jast_core_a1 () = has_array "jastrow" "jast_core_a1" ;;
let get_jastrow_jast_core_a2 () = read_float_array "jastrow" "jast_core_a2";;
let set_jastrow_jast_core_a2 = write_float_array "jastrow" "jast_core_a2";;
let has_jastrow_jast_core_a2 () = has_array "jastrow" "jast_core_a2" ;;
let get_jastrow_jast_core_b1 () = read_float_array "jastrow" "jast_core_b1";;
let set_jastrow_jast_core_b1 = write_float_array "jastrow" "jast_core_b1";;
let has_jastrow_jast_core_b1 () = has_array "jastrow" "jast_core_b1" ;;
let get_jastrow_jast_core_b2 () = read_float_array "jastrow" "jast_core_b2";;
let set_jastrow_jast_core_b2 = write_float_array "jastrow" "jast_core_b2";;
let has_jastrow_jast_core_b2 () = has_array "jastrow" "jast_core_b2" ;;
let get_jastrow_jast_1b_type () = read_int "jastrow" "jast_1b_type";;
let set_jastrow_jast_1b_type = write_int "jastrow" "jast_1b_type" ;;
let has_jastrow_jast_1b_type () = has "jastrow" "jast_1b_type" ;;
let get_jastrow_jast_1btanh_pen () = read_float_array "jastrow" "jast_1btanh_pen";;
let set_jastrow_jast_1btanh_pen = write_float_array "jastrow" "jast_1btanh_pen";;
let has_jastrow_jast_1btanh_pen () = has_array "jastrow" "jast_1btanh_pen" ;;
let get_jastrow_jast_1berf_pen () = read_float_array "jastrow" "jast_1berf_pen";;
let set_jastrow_jast_1berf_pen = write_float_array "jastrow" "jast_1berf_pen";;
let has_jastrow_jast_1berf_pen () = has_array "jastrow" "jast_1berf_pen" ;;
let get_jastrow_jast_1bgauss_pen () = read_float_array "jastrow" "jast_1bgauss_pen";;
let set_jastrow_jast_1bgauss_pen = write_float_array "jastrow" "jast_1bgauss_pen";;
let has_jastrow_jast_1bgauss_pen () = has_array "jastrow" "jast_1bgauss_pen" ;;
let get_blocks_empty () = read_int "blocks" "empty";;
let set_blocks_empty = write_int "blocks" "empty" ;;
let has_blocks_empty () = has "blocks" "empty" ;;
let get_pseudo_ao_pseudo_grid () = read_float_array "pseudo" "ao_pseudo_grid";;
let set_pseudo_ao_pseudo_grid = write_float_array "pseudo" "ao_pseudo_grid";;
let has_pseudo_ao_pseudo_grid () = has_array "pseudo" "ao_pseudo_grid" ;;
let get_pseudo_do_pseudo () = read_bool "pseudo" "do_pseudo";;
let set_pseudo_do_pseudo = write_bool "pseudo" "do_pseudo" ;;
let has_pseudo_do_pseudo () = has "pseudo" "do_pseudo" ;;
let get_pseudo_mo_pseudo_grid () = read_float_array "pseudo" "mo_pseudo_grid";;
let set_pseudo_mo_pseudo_grid = write_float_array "pseudo" "mo_pseudo_grid";;
let has_pseudo_mo_pseudo_grid () = has_array "pseudo" "mo_pseudo_grid" ;;
let get_pseudo_pseudo_dz_k () = read_float_array "pseudo" "pseudo_dz_k";;
let set_pseudo_pseudo_dz_k = write_float_array "pseudo" "pseudo_dz_k";;
let has_pseudo_pseudo_dz_k () = has_array "pseudo" "pseudo_dz_k" ;;
let get_pseudo_pseudo_dz_kl () = read_float_array "pseudo" "pseudo_dz_kl";;
let set_pseudo_pseudo_dz_kl = write_float_array "pseudo" "pseudo_dz_kl";;
let has_pseudo_pseudo_dz_kl () = has_array "pseudo" "pseudo_dz_kl" ;;
let get_pseudo_pseudo_grid_rmax () = read_float "pseudo" "pseudo_grid_rmax";;
let set_pseudo_pseudo_grid_rmax = write_float "pseudo" "pseudo_grid_rmax" ;;
let has_pseudo_pseudo_grid_rmax () = has "pseudo" "pseudo_grid_rmax" ;;
let get_pseudo_pseudo_grid_size () = read_int "pseudo" "pseudo_grid_size";;
let set_pseudo_pseudo_grid_size = write_int "pseudo" "pseudo_grid_size" ;;
let has_pseudo_pseudo_grid_size () = has "pseudo" "pseudo_grid_size" ;;
let get_pseudo_pseudo_klocmax () = read_int "pseudo" "pseudo_klocmax";;
let set_pseudo_pseudo_klocmax = write_int "pseudo" "pseudo_klocmax" ;;
let has_pseudo_pseudo_klocmax () = has "pseudo" "pseudo_klocmax" ;;
let get_pseudo_pseudo_kmax () = read_int "pseudo" "pseudo_kmax";;
let set_pseudo_pseudo_kmax = write_int "pseudo" "pseudo_kmax" ;;
let has_pseudo_pseudo_kmax () = has "pseudo" "pseudo_kmax" ;;
let get_pseudo_pseudo_lmax () = read_int "pseudo" "pseudo_lmax";;
let set_pseudo_pseudo_lmax = write_int "pseudo" "pseudo_lmax" ;;
let has_pseudo_pseudo_lmax () = has "pseudo" "pseudo_lmax" ;;
let get_pseudo_pseudo_n_k () = read_int_array "pseudo" "pseudo_n_k";;
let set_pseudo_pseudo_n_k = write_int_array "pseudo" "pseudo_n_k";;
let has_pseudo_pseudo_n_k () = has_array "pseudo" "pseudo_n_k" ;;
let get_pseudo_pseudo_n_kl () = read_int_array "pseudo" "pseudo_n_kl";;
let set_pseudo_pseudo_n_kl = write_int_array "pseudo" "pseudo_n_kl";;
let has_pseudo_pseudo_n_kl () = has_array "pseudo" "pseudo_n_kl" ;;
let get_pseudo_pseudo_v_k () = read_float_array "pseudo" "pseudo_v_k";;
let set_pseudo_pseudo_v_k = write_float_array "pseudo" "pseudo_v_k";;
let has_pseudo_pseudo_v_k () = has_array "pseudo" "pseudo_v_k" ;;
let get_pseudo_pseudo_v_kl () = read_float_array "pseudo" "pseudo_v_kl";;
let set_pseudo_pseudo_v_kl = write_float_array "pseudo" "pseudo_v_kl";;
let has_pseudo_pseudo_v_kl () = has_array "pseudo" "pseudo_v_kl" ;;
let get_properties_ci_dress () = read_bool "properties" "ci_dress";;
let set_properties_ci_dress = write_bool "properties" "ci_dress" ;;
let has_properties_ci_dress () = has "properties" "ci_dress" ;;
let get_properties_ci_dress_mu () = read_bool "properties" "ci_dress_mu";;
let set_properties_ci_dress_mu = write_bool "properties" "ci_dress_mu" ;;
let has_properties_ci_dress_mu () = has "properties" "ci_dress_mu" ;;
let get_properties_ci_h_matrix () = read_bool "properties" "ci_h_matrix";;
let set_properties_ci_h_matrix = write_bool "properties" "ci_h_matrix" ;;
let has_properties_ci_h_matrix () = has "properties" "ci_h_matrix" ;;
let get_properties_ci_h_matrix_diag () = read_bool "properties" "ci_h_matrix_diag";;
let set_properties_ci_h_matrix_diag = write_bool "properties" "ci_h_matrix_diag" ;;
let has_properties_ci_h_matrix_diag () = has "properties" "ci_h_matrix_diag" ;;
let get_properties_ci_h_psidet () = read_bool "properties" "ci_h_psidet";;
let set_properties_ci_h_psidet = write_bool "properties" "ci_h_psidet" ;;
let has_properties_ci_h_psidet () = has "properties" "ci_h_psidet" ;;
let get_properties_ci_h_transcor_psi () = read_bool "properties" "ci_h_transcor_psi";;
let set_properties_ci_h_transcor_psi = write_bool "properties" "ci_h_transcor_psi" ;;
let has_properties_ci_h_transcor_psi () = has "properties" "ci_h_transcor_psi" ;;
let get_properties_ci_overlap_matrix () = read_bool "properties" "ci_overlap_matrix";;
let set_properties_ci_overlap_matrix = write_bool "properties" "ci_overlap_matrix" ;;
let has_properties_ci_overlap_matrix () = has "properties" "ci_overlap_matrix" ;;
let get_properties_ci_overlap_psidet () = read_bool "properties" "ci_overlap_psidet";;
let set_properties_ci_overlap_psidet = write_bool "properties" "ci_overlap_psidet" ;;
let has_properties_ci_overlap_psidet () = has "properties" "ci_overlap_psidet" ;;
let get_properties_dipole () = read_bool "properties" "dipole";;
let set_properties_dipole = write_bool "properties" "dipole" ;;
let has_properties_dipole () = has "properties" "dipole" ;;
let get_properties_drift_mod () = read_bool "properties" "drift_mod";;
let set_properties_drift_mod = write_bool "properties" "drift_mod" ;;
let has_properties_drift_mod () = has "properties" "drift_mod" ;;
let get_properties_e_kin () = read_bool "properties" "e_kin";;
let set_properties_e_kin = write_bool "properties" "e_kin" ;;
let has_properties_e_kin () = has "properties" "e_kin" ;;
let get_properties_e_loc () = read_bool "properties" "e_loc";;
let set_properties_e_loc = write_bool "properties" "e_loc" ;;
let has_properties_e_loc () = has "properties" "e_loc" ;;
let get_properties_e_loc_zv () = read_bool "properties" "e_loc_zv";;
let set_properties_e_loc_zv = write_bool "properties" "e_loc_zv" ;;
let has_properties_e_loc_zv () = has "properties" "e_loc_zv" ;;
let get_properties_e_nucl () = read_bool "properties" "e_nucl";;
let set_properties_e_nucl = write_bool "properties" "e_nucl" ;;
let has_properties_e_nucl () = has "properties" "e_nucl" ;;
let get_properties_e_nucl_elec () = read_bool "properties" "e_nucl_elec";;
let set_properties_e_nucl_elec = write_bool "properties" "e_nucl_elec" ;;
let has_properties_e_nucl_elec () = has "properties" "e_nucl_elec" ;;
let get_properties_e_pot () = read_bool "properties" "e_pot";;
let set_properties_e_pot = write_bool "properties" "e_pot" ;;
let has_properties_e_pot () = has "properties" "e_pot" ;;
let get_properties_eff_pot_deriv_mu () = read_bool "properties" "eff_pot_deriv_mu";;
let set_properties_eff_pot_deriv_mu = write_bool "properties" "eff_pot_deriv_mu" ;;
let has_properties_eff_pot_deriv_mu () = has "properties" "eff_pot_deriv_mu" ;;
let get_properties_eff_pot_deriv_mu_elec () = read_bool "properties" "eff_pot_deriv_mu_elec";;
let set_properties_eff_pot_deriv_mu_elec = write_bool "properties" "eff_pot_deriv_mu_elec" ;;
let has_properties_eff_pot_deriv_mu_elec () = has "properties" "eff_pot_deriv_mu_elec" ;;
let get_properties_eff_pot_mu () = read_bool "properties" "eff_pot_mu";;
let set_properties_eff_pot_mu = write_bool "properties" "eff_pot_mu" ;;
let has_properties_eff_pot_mu () = has "properties" "eff_pot_mu" ;;
let get_properties_eff_pot_mu_elec () = read_bool "properties" "eff_pot_mu_elec";;
let set_properties_eff_pot_mu_elec = write_bool "properties" "eff_pot_mu_elec" ;;
let has_properties_eff_pot_mu_elec () = has "properties" "eff_pot_mu_elec" ;;
let get_properties_eff_pot_mu_simple () = read_bool "properties" "eff_pot_mu_simple";;
let set_properties_eff_pot_mu_simple = write_bool "properties" "eff_pot_mu_simple" ;;
let has_properties_eff_pot_mu_simple () = has "properties" "eff_pot_mu_simple" ;;
let get_properties_emudiff () = read_bool "properties" "emudiff";;
let set_properties_emudiff = write_bool "properties" "emudiff" ;;
let has_properties_emudiff () = has "properties" "emudiff" ;;
let get_properties_energy_mu () = read_bool "properties" "energy_mu";;
let set_properties_energy_mu = write_bool "properties" "energy_mu" ;;
let has_properties_energy_mu () = has "properties" "energy_mu" ;;
let get_properties_pop_weight () = read_bool "properties" "pop_weight";;
let set_properties_pop_weight = write_bool "properties" "pop_weight" ;;
let has_properties_pop_weight () = has "properties" "pop_weight" ;;
let get_properties_psi_norm () = read_bool "properties" "psi_norm";;
let set_properties_psi_norm = write_bool "properties" "psi_norm" ;;
let has_properties_psi_norm () = has "properties" "psi_norm" ;;
let get_properties_three_body_mu () = read_bool "properties" "three_body_mu";;
let set_properties_three_body_mu = write_bool "properties" "three_body_mu" ;;
let has_properties_three_body_mu () = has "properties" "three_body_mu" ;;
let get_properties_wf_extension () = read_bool "properties" "wf_extension";;
let set_properties_wf_extension = write_bool "properties" "wf_extension" ;;
let has_properties_wf_extension () = has "properties" "wf_extension" ;;
let get_ezfio_creation () = read_string "ezfio" "creation";;
let set_ezfio_creation = write_string "ezfio" "creation" ;;
let has_ezfio_creation () = has "ezfio" "creation" ;;
let get_ezfio_user () = read_string "ezfio" "user";;
let set_ezfio_user = write_string "ezfio" "user" ;;
let has_ezfio_user () = has "ezfio" "user" ;;
let get_ezfio_library () = read_string "ezfio" "library";;
let set_ezfio_library = write_string "ezfio" "library" ;;
let has_ezfio_library () = has "ezfio" "library" ;;
let get_ezfio_last_library () = read_string "ezfio" "last_library";;
let set_ezfio_last_library = write_string "ezfio" "last_library" ;;
let has_ezfio_last_library () = has "ezfio" "last_library" ;;
