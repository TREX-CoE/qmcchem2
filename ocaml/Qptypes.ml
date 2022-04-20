
let bytes_of_int64 i =
  let result = Bytes.create 8 in
  Bytes.set_int64_ne result 0 i;
  result

let bytes_of_int i =
  Int64.of_int i
  |> bytes_of_int64


let int64_of_bytes b =
  Bytes.get_int64_ne b 0


let int_of_bytes b =
  int64_of_bytes b
  |> Int64.to_int


let float_of_bytes b =
  int64_of_bytes b
  |> Int64.float_of_bits


let bytes_of_float f =
  Int64.bits_of_float f
  |> bytes_of_int64


open Sexplib.Std
let warning = print_string

module Compute_node : sig
  type t [@@deriving sexp]
  val to_string : t -> string
  val of_string :  string -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = string [@@deriving sexp]
  let to_string x = x
  let of_string  x = ( assert (x <> "") ; x )
  let to_string x =  x
  let to_bytes  x = Bytes.of_string  x
  let of_bytes  b = Bytes.to_string  b
end


module Block_id : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x > 0) ; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Weight : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( assert (x >= 0.) ; x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Rst_string : sig
  type t [@@deriving sexp]
  val to_string : t -> string
  val of_string :  string -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = string [@@deriving sexp]
  let to_string x = x
  let of_string  x = (  x )
  let to_string x =  x
  let to_bytes  x = Bytes.of_string  x
  let of_bytes  b = Bytes.to_string  b
end


module MD5 : sig
  type t [@@deriving sexp]
  val to_string : t -> string
  val of_string :  string -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = string [@@deriving sexp]
  let to_string x = x
  let of_string  x = ( assert ((String.length x) = 32);
  assert (
    let a =
      Array.init (String.length x) (fun i -> x.[i])
    in
    Array.fold_left (fun accu x -> accu && (x < 'g')) true a
    ); x )
  let to_string x =  x
  let to_bytes  x = Bytes.of_string  x
  let of_bytes  b = Bytes.to_string  b
end


module Elec_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x > 0) ; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Elec_beta_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x >= 0) ; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Elec_alpha_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x > 0) ; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module PT2_energy : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( assert (x >=0.) ; x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Threshold : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( assert (x >= 0.) ;
  assert (x <= 1.) ; x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module AO_prim_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x > 0) ; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module AO_expo : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if (x < 0.) then
    raise (Invalid_argument (Printf.sprintf "AO_expo : (x >= 0.) : x=%f" x)); x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module AO_coef : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = (  x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module MO_occ : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if x < 0. then 0.  else
  if x > 2. then 2.  else x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module MO_coef : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = (  x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Bitmask_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x > 0) ; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Bit_kind : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( begin match x with
  | 1 | 2 | 4 | 8 -> ()
  | _ -> raise (Invalid_argument "Bit_kind should be (1|2|4|8).")
  end; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Bit_kind_size : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( begin match x with
  | 8 | 16 | 32 | 64 -> ()
  | _ -> raise (Invalid_argument "Bit_kind_size should be (8|16|32|64).")
  end; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module States_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x > 0) ;
  if (x > 1000) then
    warning "More than 1000 states"; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Det_number_max : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( assert (x > 0) ;
  if (x > 100_000_000) then
    warning "More than 100 million determinants"; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Non_empty_string : sig
  type t [@@deriving sexp]
  val to_string : t -> string
  val of_string :  string -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = string [@@deriving sexp]
  let to_string x = x
  let of_string  x = ( if (x = "") then
    raise (Invalid_argument "Non_empty_string"); x )
  let to_string x =  x
  let to_bytes  x = Bytes.of_string  x
  let of_bytes  b = Bytes.to_string  b
end


module Strictly_negative_int : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( if not (x < 0) then
    raise (Invalid_argument (Printf.sprintf "Strictly_negative_int : (x < 0) : x=%d" x)); x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Normalized_float : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if (x < 0.) || (x > 1.) then
    raise (Invalid_argument (Printf.sprintf "Normalized_float : (0. <= x <= 1.) : x=%f" x)); x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Det_coef : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if (x < -1.) || (x > 1.) then
    raise (Invalid_argument (Printf.sprintf "Det_coef : (-1. <= x <= 1.) : x=%f" x)); x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Negative_int : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( if not (x <= 0) then
    raise (Invalid_argument (Printf.sprintf "Negative_int : (x <= 0) : x=%d" x));
  assert (x <= 0) ; x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Strictly_positive_int : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( if not (x > 0) then
    raise (Invalid_argument (Printf.sprintf "Strictly_positive_int : (x > 0) : x=%d" x)); x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Positive_int : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val of_int :  int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int [@@deriving sexp]
  let to_int x = x
  let of_int  x = ( if not (x >= 0) then
    raise (Invalid_argument (Printf.sprintf "Positive_int : (x >= 0) : x=%d" x)); x )
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let of_bytes  b = int_of_bytes  b
end


module Positive_int64 : sig
  type t [@@deriving sexp]
  val to_int64 : t -> int64
  val of_int64 :  int64 -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = int64 [@@deriving sexp]
  let to_int64 x = x
  let of_int64  x = ( if not (x >= 0L) then
    raise (Invalid_argument (Printf.sprintf "Positive_int64 : (x >= 0L) : x=%s" (Int64.to_string x))); x )
  let to_string x = Int64.to_string x
  let to_bytes  x = bytes_of_int64  x
  let of_bytes  b = int64_of_bytes  b
end


module Strictly_negative_float : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if not (x < 0.) then
    raise (Invalid_argument (Printf.sprintf "Strictly_negative_float : (x < 0.) : x=%f" x)); x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Negative_float : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if not (x <= 0.) then
    raise (Invalid_argument (Printf.sprintf "Negative_float : (x <= 0.) : x=%f" x)); x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Strictly_positive_float : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if not (x > 0.) then
    raise (Invalid_argument (Printf.sprintf "Strictly_positive_float : (x > 0.) : x=%f" x)); x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module Positive_float : sig
  type t [@@deriving sexp]
  val to_float : t -> float
  val of_float :  float -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = float [@@deriving sexp]
  let to_float x = x
  let of_float  x = ( if not (x >= 0.) then
    raise (Invalid_argument (Printf.sprintf "Positive_float : (x >= 0.) : x=%f"  x)); x )
  let to_string x = string_of_float x
  let to_bytes  x = bytes_of_float  x
  let of_bytes  b = float_of_bytes  b
end


module MO_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val get_max : unit -> int
  val of_int : ?min:int -> ?max:int -> int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
end = struct
  type t = int [@@deriving sexp]
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let get_max () =
    if (Ezfio.has_mo_basis_mo_num ()) then
      Ezfio.get_mo_basis_mo_num ()
    else
      10_000
  let get_min () =
      1
  let to_int x = x
  let of_int ?(min=get_min ()) ?(max=get_max ()) x =
    begin
      assert (x >= min) ;
      if (x > 10_000) then
        warning "More than 10_000 MOs";
      begin
        match max with
        | 1 -> ()
        | i  ->
          if ( x > i ) then
            raise (Invalid_argument (Printf.sprintf "MO_number: %s" (string_of_int x) ))
      end ;
      x
    end
end

module AO_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val get_max : unit -> int
  val of_int : ?min:int -> ?max:int -> int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
end = struct
  type t = int [@@deriving sexp]
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let get_max () =
    if (Ezfio.has_ao_basis_ao_num ()) then
      Ezfio.get_ao_basis_ao_num ()
    else
      10_000
  let get_min () =
      1
  let to_int x = x
  let of_int ?(min=get_min ()) ?(max=get_max ()) x =
    begin
      assert (x >= min) ;
      if (x > 10_000) then
        warning "More than 10_000 AOs";
      begin
        match max with
        | 1 -> ()
        | i  ->
          if ( x > i ) then
            raise (Invalid_argument (Printf.sprintf "AO_number: %s" (string_of_int x) ))
      end ;
      x
    end
end

module Nucl_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val get_max : unit -> int
  val of_int : ?min:int -> ?max:int -> int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
end = struct
  type t = int [@@deriving sexp]
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let get_max () =
    if (Ezfio.has_nuclei_nucl_num ()) then
      Ezfio.get_nuclei_nucl_num ()
    else
      10_000
  let get_min () =
      1
  let to_int x = x
  let of_int ?(min=get_min ()) ?(max=get_max ()) x =
    begin
      assert (x >= min) ;
      if (x > 10_000) then
        warning "More than 10_000 nuclei";
      begin
        match max with
        | 1 -> ()
        | i  ->
          if ( x > i ) then
            raise (Invalid_argument (Printf.sprintf "Nucl_number: %s" (string_of_int x) ))
      end ;
      x
    end
end

module N_int_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val get_max : unit -> int
  val of_int : ?min:int -> ?max:int -> int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
end = struct
  type t = int [@@deriving sexp]
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let get_max () =
    if (Ezfio.has_spindeterminants_n_int ()) then
      Ezfio.get_spindeterminants_n_int ()
    else
      30
  let get_min () =
      1
  let to_int x = x
  let of_int ?(min=get_min ()) ?(max=get_max ()) x =
    begin
      assert (x >= min) ;
      if (x > 30) then
        warning "N_int > 30";
      begin
        match max with
        | 1 -> ()
        | i  ->
          if ( x > i ) then
            raise (Invalid_argument (Printf.sprintf "N_int_number: %s" (string_of_int x) ))
      end ;
      x
    end
end

module Det_number : sig
  type t [@@deriving sexp]
  val to_int : t -> int
  val get_max : unit -> int
  val of_int : ?min:int -> ?max:int -> int -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
end = struct
  type t = int [@@deriving sexp]
  let to_string x = string_of_int x
  let to_bytes  x = bytes_of_int  x
  let get_max () =
    if (Ezfio.has_spindeterminants_n_det ()) then
      Ezfio.get_spindeterminants_n_det ()
    else
      100_000_000
  let get_min () =
      1
  let to_int x = x
  let of_int ?(min=get_min ()) ?(max=get_max ()) x =
    begin
      assert (x >= min) ;
      if (x > 100_000_000) then
        warning "More than 100 million determinants";
      begin
        match max with
        | 1 -> ()
        | i  ->
          if ( x > i ) then
            raise (Invalid_argument (Printf.sprintf "Det_number: %s" (string_of_int x) ))
      end ;
      x
    end
end
let decode_ezfio_message msg =
match msg with 
 | "get_ezfio_user" -> Ezfio.read_string "ezfio" "user"
 | "get_ezfio_library" -> Ezfio.read_string "ezfio" "library"
 | "get_ezfio_last_library" -> Ezfio.read_string "ezfio" "last_library"
 | "get_properties_e_nucl" -> Ezfio.read_string "properties" "e_nucl"
 | "get_properties_e_pot" -> Ezfio.read_string "properties" "e_pot"
 | "get_properties_e_kin" -> Ezfio.read_string "properties" "e_kin"
 | "get_properties_e_loc" -> Ezfio.read_string "properties" "e_loc"
 | "get_properties_e_loc_zv" -> Ezfio.read_string "properties" "e_loc_zv"
 | "get_properties_dipole" -> Ezfio.read_string "properties" "dipole"
 | "get_properties_wf_extension" -> Ezfio.read_string "properties" "wf_extension"
 | "get_properties_pop_weight" -> Ezfio.read_string "properties" "pop_weight"
 | "get_properties_drift_mod" -> Ezfio.read_string "properties" "drift_mod"
 | "get_properties_ci_h_transcor_psi" -> Ezfio.read_string "properties" "ci_h_transcor_psi"
 | "get_properties_emudiff" -> Ezfio.read_string "properties" "emudiff"
 | "get_properties_psi_norm" -> Ezfio.read_string "properties" "psi_norm"
 | "get_properties_ci_dress" -> Ezfio.read_string "properties" "ci_dress"
 | "get_properties_ci_dress_mu" -> Ezfio.read_string "properties" "ci_dress_mu"
 | "get_properties_ci_overlap_psidet" -> Ezfio.read_string "properties" "ci_overlap_psidet"
 | "get_properties_ci_h_psidet" -> Ezfio.read_string "properties" "ci_h_psidet"
 | "get_properties_ci_overlap_matrix" -> Ezfio.read_string "properties" "ci_overlap_matrix"
 | "get_properties_ci_h_matrix" -> Ezfio.read_string "properties" "ci_h_matrix"
 | "get_properties_ci_h_matrix_diag" -> Ezfio.read_string "properties" "ci_h_matrix_diag"
 | "get_properties_ci_dress_mu_opt" -> Ezfio.read_string "properties" "ci_dress_mu_opt"
 | "get_ao_basis_ao_num" -> Ezfio.read_string "ao_basis" "ao_num"
 | "get_ao_basis_ao_prim_num" ->
             Ezfio.read_string_array "ao_basis" "ao_prim_num"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_ao_basis_ao_nucl" ->
             Ezfio.read_string_array "ao_basis" "ao_nucl"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_ao_basis_ao_power" ->
             Ezfio.read_string_array "ao_basis" "ao_power"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_ao_basis_ao_coef" ->
             Ezfio.read_string_array "ao_basis" "ao_coef"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_ao_basis_ao_expo" ->
             Ezfio.read_string_array "ao_basis" "ao_expo"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_mo_basis_mo_num" -> Ezfio.read_string "mo_basis" "mo_num"
 | "get_mo_basis_mo_coef" ->
             Ezfio.read_string_array "mo_basis" "mo_coef"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_mo_basis_mo_classif" ->
             Ezfio.read_string_array "mo_basis" "mo_classif"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_mo_basis_mo_energy" ->
             Ezfio.read_string_array "mo_basis" "mo_energy"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_mo_basis_mo_occ" ->
             Ezfio.read_string_array "mo_basis" "mo_occ"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_mo_basis_mo_symmetry" ->
             Ezfio.read_string_array "mo_basis" "mo_symmetry"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_electrons_elec_alpha_num" -> Ezfio.read_string "electrons" "elec_alpha_num"
 | "get_electrons_elec_beta_num" -> Ezfio.read_string "electrons" "elec_beta_num"
 | "get_electrons_elec_walk_num_tot" -> Ezfio.read_string "electrons" "elec_walk_num_tot"
 | "get_electrons_elec_walk_num" -> Ezfio.read_string "electrons" "elec_walk_num"
 | "get_electrons_elec_coord_pool" ->
             Ezfio.read_string_array "electrons" "elec_coord_pool"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_electrons_elec_coord_pool_size" -> Ezfio.read_string "electrons" "elec_coord_pool_size"
 | "get_electrons_elec_fitcusp_radius" -> Ezfio.read_string "electrons" "elec_fitcusp_radius"
 | "get_nuclei_nucl_num" -> Ezfio.read_string "nuclei" "nucl_num"
 | "get_nuclei_nucl_label" ->
             Ezfio.read_string_array "nuclei" "nucl_label"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_nuclei_nucl_charge" ->
             Ezfio.read_string_array "nuclei" "nucl_charge"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_nuclei_nucl_coord" ->
             Ezfio.read_string_array "nuclei" "nucl_coord"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_n_det_alpha" -> Ezfio.read_string "spindeterminants" "n_det_alpha"
 | "get_spindeterminants_n_det_beta" -> Ezfio.read_string "spindeterminants" "n_det_beta"
 | "get_spindeterminants_n_det" -> Ezfio.read_string "spindeterminants" "n_det"
 | "get_spindeterminants_n_int" -> Ezfio.read_string "spindeterminants" "n_int"
 | "get_spindeterminants_bit_kind" -> Ezfio.read_string "spindeterminants" "bit_kind"
 | "get_spindeterminants_n_states" -> Ezfio.read_string "spindeterminants" "n_states"
 | "get_spindeterminants_psi_det_alpha" ->
             Ezfio.read_string_array "spindeterminants" "psi_det_alpha"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_det_beta" ->
             Ezfio.read_string_array "spindeterminants" "psi_det_beta"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_coef_matrix_rows" ->
             Ezfio.read_string_array "spindeterminants" "psi_coef_matrix_rows"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_coef_matrix_columns" ->
             Ezfio.read_string_array "spindeterminants" "psi_coef_matrix_columns"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_coef_matrix_values" ->
             Ezfio.read_string_array "spindeterminants" "psi_coef_matrix_values"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_n_svd_coefs_unique" -> Ezfio.read_string "spindeterminants" "n_svd_coefs_unique"
 | "get_spindeterminants_n_svd_coefs" -> Ezfio.read_string "spindeterminants" "n_svd_coefs"
 | "get_spindeterminants_n_svd_selected" -> Ezfio.read_string "spindeterminants" "n_svd_selected"
 | "get_spindeterminants_n_svd_toselect" -> Ezfio.read_string "spindeterminants" "n_svd_toselect"
 | "get_spindeterminants_psi_svd_alpha_unique" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_alpha_unique"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_beta_unique" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_beta_unique"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_coefs_unique" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_coefs_unique"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_alpha" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_alpha"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_beta" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_beta"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_coefs" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_coefs"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_alpha_numselected" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_alpha_numselected"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_beta_numselected" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_beta_numselected"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_alpha_numtoselect" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_alpha_numtoselect"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_spindeterminants_psi_svd_beta_numtoselect" ->
             Ezfio.read_string_array "spindeterminants" "psi_svd_beta_numtoselect"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_simulation_do_run" -> Ezfio.read_string "simulation" "do_run"
 | "get_simulation_stop_time" -> Ezfio.read_string "simulation" "stop_time"
 | "get_simulation_equilibration" -> Ezfio.read_string "simulation" "equilibration"
 | "get_simulation_http_server" -> Ezfio.read_string "simulation" "http_server"
 | "get_simulation_do_jast" -> Ezfio.read_string "simulation" "do_jast"
 | "get_simulation_nucl_fitcusp_factor" -> Ezfio.read_string "simulation" "nucl_fitcusp_factor"
 | "get_simulation_method" -> Ezfio.read_string "simulation" "method"
 | "get_simulation_block_time" -> Ezfio.read_string "simulation" "block_time"
 | "get_simulation_sampling" -> Ezfio.read_string "simulation" "sampling"
 | "get_simulation_save_data" -> Ezfio.read_string "simulation" "save_data"
 | "get_simulation_time_step" -> Ezfio.read_string "simulation" "time_step"
 | "get_simulation_print_level" -> Ezfio.read_string "simulation" "print_level"
 | "get_simulation_ci_threshold" -> Ezfio.read_string "simulation" "ci_threshold"
 | "get_simulation_md5_key" -> Ezfio.read_string "simulation" "md5_key"
 | "get_simulation_e_ref" -> Ezfio.read_string "simulation" "e_ref"
 | "get_simulation_e_trial" -> Ezfio.read_string "simulation" "e_trial"
 | "get_simulation_srmc_projection_time" -> Ezfio.read_string "simulation" "srmc_projection_time"
 | "get_simulation_use_trexio" -> Ezfio.read_string "simulation" "use_trexio"
 | "get_simulation_use_qmckl" -> Ezfio.read_string "simulation" "use_qmckl"
 | "get_simulation_trexio_filename" -> Ezfio.read_string "simulation" "trexio_filename"
 | "get_jastrow_jast_type" -> Ezfio.read_string "jastrow" "jast_type"
 | "get_jastrow_jast_a_up_up" -> Ezfio.read_string "jastrow" "jast_a_up_up"
 | "get_jastrow_jast_a_up_dn" -> Ezfio.read_string "jastrow" "jast_a_up_dn"
 | "get_jastrow_jast_b_up_up" -> Ezfio.read_string "jastrow" "jast_b_up_up"
 | "get_jastrow_jast_b_up_dn" -> Ezfio.read_string "jastrow" "jast_b_up_dn"
 | "get_jastrow_mu_erf" -> Ezfio.read_string "jastrow" "mu_erf"
 | "get_jastrow_jast_pen" ->
             Ezfio.read_string_array "jastrow" "jast_pen"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_een_e_a" ->
             Ezfio.read_string_array "jastrow" "jast_een_e_a"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_een_e_b" ->
             Ezfio.read_string_array "jastrow" "jast_een_e_b"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_een_n" ->
             Ezfio.read_string_array "jastrow" "jast_een_n"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_core_a1" ->
             Ezfio.read_string_array "jastrow" "jast_core_a1"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_core_a2" ->
             Ezfio.read_string_array "jastrow" "jast_core_a2"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_core_b1" ->
             Ezfio.read_string_array "jastrow" "jast_core_b1"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_core_b2" ->
             Ezfio.read_string_array "jastrow" "jast_core_b2"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_1b_type" -> Ezfio.read_string "jastrow" "jast_1b_type"
 | "get_jastrow_jast_1btanh_pen" ->
             Ezfio.read_string_array "jastrow" "jast_1btanh_pen"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_1berf_pen" ->
             Ezfio.read_string_array "jastrow" "jast_1berf_pen"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_jastrow_jast_1bgauss_pen" ->
             Ezfio.read_string_array "jastrow" "jast_1bgauss_pen"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_blocks_empty" -> Ezfio.read_string "blocks" "empty"
 | "get_pseudo_ao_pseudo_grid" ->
             Ezfio.read_string_array "pseudo" "ao_pseudo_grid"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_pseudo_do_pseudo" -> Ezfio.read_string "pseudo" "do_pseudo"
 | "get_pseudo_mo_pseudo_grid" ->
             Ezfio.read_string_array "pseudo" "mo_pseudo_grid"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_pseudo_pseudo_dz_k" ->
             Ezfio.read_string_array "pseudo" "pseudo_dz_k"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_pseudo_pseudo_dz_kl" ->
             Ezfio.read_string_array "pseudo" "pseudo_dz_kl"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_pseudo_pseudo_grid_rmax" -> Ezfio.read_string "pseudo" "pseudo_grid_rmax"
 | "get_pseudo_pseudo_grid_size" -> Ezfio.read_string "pseudo" "pseudo_grid_size"
 | "get_pseudo_pseudo_klocmax" -> Ezfio.read_string "pseudo" "pseudo_klocmax"
 | "get_pseudo_pseudo_kmax" -> Ezfio.read_string "pseudo" "pseudo_kmax"
 | "get_pseudo_pseudo_lmax" -> Ezfio.read_string "pseudo" "pseudo_lmax"
 | "get_pseudo_pseudo_n_k" ->
             Ezfio.read_string_array "pseudo" "pseudo_n_k"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_pseudo_pseudo_n_kl" ->
             Ezfio.read_string_array "pseudo" "pseudo_n_kl"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_pseudo_pseudo_v_k" ->
             Ezfio.read_string_array "pseudo" "pseudo_v_k"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "get_pseudo_pseudo_v_kl" ->
             Ezfio.read_string_array "pseudo" "pseudo_v_kl"
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat " "
 | "has_ezfio_creation" -> if (Ezfio.has_ezfio_creation ()) then "T" else "F"
 | "has_ezfio_user" -> if (Ezfio.has_ezfio_user ()) then "T" else "F"
 | "has_ezfio_library" -> if (Ezfio.has_ezfio_library ()) then "T" else "F"
 | "has_ezfio_last_library" -> if (Ezfio.has_ezfio_last_library ()) then "T" else "F"
 | "has_properties_e_nucl" -> if (Ezfio.has_properties_e_nucl ()) then "T" else "F"
 | "has_properties_e_pot" -> if (Ezfio.has_properties_e_pot ()) then "T" else "F"
 | "has_properties_e_kin" -> if (Ezfio.has_properties_e_kin ()) then "T" else "F"
 | "has_properties_e_loc" -> if (Ezfio.has_properties_e_loc ()) then "T" else "F"
 | "has_properties_e_loc_zv" -> if (Ezfio.has_properties_e_loc_zv ()) then "T" else "F"
 | "has_properties_dipole" -> if (Ezfio.has_properties_dipole ()) then "T" else "F"
 | "has_properties_wf_extension" -> if (Ezfio.has_properties_wf_extension ()) then "T" else "F"
 | "has_properties_pop_weight" -> if (Ezfio.has_properties_pop_weight ()) then "T" else "F"
 | "has_properties_drift_mod" -> if (Ezfio.has_properties_drift_mod ()) then "T" else "F"
 | "has_properties_ci_h_transcor_psi" -> if (Ezfio.has_properties_ci_h_transcor_psi ()) then "T" else "F"
 | "has_properties_emudiff" -> if (Ezfio.has_properties_emudiff ()) then "T" else "F"
 | "has_properties_psi_norm" -> if (Ezfio.has_properties_psi_norm ()) then "T" else "F"
 | "has_properties_ci_dress" -> if (Ezfio.has_properties_ci_dress ()) then "T" else "F"
 | "has_properties_ci_dress_mu" -> if (Ezfio.has_properties_ci_dress_mu ()) then "T" else "F"
 | "has_properties_ci_overlap_psidet" -> if (Ezfio.has_properties_ci_overlap_psidet ()) then "T" else "F"
 | "has_properties_ci_h_psidet" -> if (Ezfio.has_properties_ci_h_psidet ()) then "T" else "F"
 | "has_properties_ci_overlap_matrix" -> if (Ezfio.has_properties_ci_overlap_matrix ()) then "T" else "F"
 | "has_properties_ci_h_matrix" -> if (Ezfio.has_properties_ci_h_matrix ()) then "T" else "F"
 | "has_properties_ci_h_matrix_diag" -> if (Ezfio.has_properties_ci_h_matrix_diag ()) then "T" else "F"
 | "has_properties_ci_dress_mu_opt" -> if (Ezfio.has_properties_ci_dress_mu_opt ()) then "T" else "F"
 | "has_ao_basis_ao_num" -> if (Ezfio.has_ao_basis_ao_num ()) then "T" else "F"
 | "has_ao_basis_ao_prim_num" -> if (Ezfio.has_ao_basis_ao_prim_num ()) then "T" else "F"
 | "has_ao_basis_ao_nucl" -> if (Ezfio.has_ao_basis_ao_nucl ()) then "T" else "F"
 | "has_ao_basis_ao_power" -> if (Ezfio.has_ao_basis_ao_power ()) then "T" else "F"
 | "has_ao_basis_ao_coef" -> if (Ezfio.has_ao_basis_ao_coef ()) then "T" else "F"
 | "has_ao_basis_ao_expo" -> if (Ezfio.has_ao_basis_ao_expo ()) then "T" else "F"
 | "has_mo_basis_mo_num" -> if (Ezfio.has_mo_basis_mo_num ()) then "T" else "F"
 | "has_mo_basis_mo_coef" -> if (Ezfio.has_mo_basis_mo_coef ()) then "T" else "F"
 | "has_mo_basis_mo_classif" -> if (Ezfio.has_mo_basis_mo_classif ()) then "T" else "F"
 | "has_mo_basis_mo_energy" -> if (Ezfio.has_mo_basis_mo_energy ()) then "T" else "F"
 | "has_mo_basis_mo_occ" -> if (Ezfio.has_mo_basis_mo_occ ()) then "T" else "F"
 | "has_mo_basis_mo_symmetry" -> if (Ezfio.has_mo_basis_mo_symmetry ()) then "T" else "F"
 | "has_electrons_elec_alpha_num" -> if (Ezfio.has_electrons_elec_alpha_num ()) then "T" else "F"
 | "has_electrons_elec_beta_num" -> if (Ezfio.has_electrons_elec_beta_num ()) then "T" else "F"
 | "has_electrons_elec_walk_num_tot" -> if (Ezfio.has_electrons_elec_walk_num_tot ()) then "T" else "F"
 | "has_electrons_elec_walk_num" -> if (Ezfio.has_electrons_elec_walk_num ()) then "T" else "F"
 | "has_electrons_elec_coord_pool" -> if (Ezfio.has_electrons_elec_coord_pool ()) then "T" else "F"
 | "has_electrons_elec_coord_pool_size" -> if (Ezfio.has_electrons_elec_coord_pool_size ()) then "T" else "F"
 | "has_electrons_elec_fitcusp_radius" -> if (Ezfio.has_electrons_elec_fitcusp_radius ()) then "T" else "F"
 | "has_nuclei_nucl_num" -> if (Ezfio.has_nuclei_nucl_num ()) then "T" else "F"
 | "has_nuclei_nucl_label" -> if (Ezfio.has_nuclei_nucl_label ()) then "T" else "F"
 | "has_nuclei_nucl_charge" -> if (Ezfio.has_nuclei_nucl_charge ()) then "T" else "F"
 | "has_nuclei_nucl_coord" -> if (Ezfio.has_nuclei_nucl_coord ()) then "T" else "F"
 | "has_spindeterminants_n_det_alpha" -> if (Ezfio.has_spindeterminants_n_det_alpha ()) then "T" else "F"
 | "has_spindeterminants_n_det_beta" -> if (Ezfio.has_spindeterminants_n_det_beta ()) then "T" else "F"
 | "has_spindeterminants_n_det" -> if (Ezfio.has_spindeterminants_n_det ()) then "T" else "F"
 | "has_spindeterminants_n_int" -> if (Ezfio.has_spindeterminants_n_int ()) then "T" else "F"
 | "has_spindeterminants_bit_kind" -> if (Ezfio.has_spindeterminants_bit_kind ()) then "T" else "F"
 | "has_spindeterminants_n_states" -> if (Ezfio.has_spindeterminants_n_states ()) then "T" else "F"
 | "has_spindeterminants_psi_det_alpha" -> if (Ezfio.has_spindeterminants_psi_det_alpha ()) then "T" else "F"
 | "has_spindeterminants_psi_det_beta" -> if (Ezfio.has_spindeterminants_psi_det_beta ()) then "T" else "F"
 | "has_spindeterminants_psi_coef_matrix_rows" -> if (Ezfio.has_spindeterminants_psi_coef_matrix_rows ()) then "T" else "F"
 | "has_spindeterminants_psi_coef_matrix_columns" -> if (Ezfio.has_spindeterminants_psi_coef_matrix_columns ()) then "T" else "F"
 | "has_spindeterminants_psi_coef_matrix_values" -> if (Ezfio.has_spindeterminants_psi_coef_matrix_values ()) then "T" else "F"
 | "has_spindeterminants_n_svd_coefs_unique" -> if (Ezfio.has_spindeterminants_n_svd_coefs_unique ()) then "T" else "F"
 | "has_spindeterminants_n_svd_coefs" -> if (Ezfio.has_spindeterminants_n_svd_coefs ()) then "T" else "F"
 | "has_spindeterminants_n_svd_selected" -> if (Ezfio.has_spindeterminants_n_svd_selected ()) then "T" else "F"
 | "has_spindeterminants_n_svd_toselect" -> if (Ezfio.has_spindeterminants_n_svd_toselect ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_alpha_unique" -> if (Ezfio.has_spindeterminants_psi_svd_alpha_unique ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_beta_unique" -> if (Ezfio.has_spindeterminants_psi_svd_beta_unique ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_coefs_unique" -> if (Ezfio.has_spindeterminants_psi_svd_coefs_unique ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_alpha" -> if (Ezfio.has_spindeterminants_psi_svd_alpha ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_beta" -> if (Ezfio.has_spindeterminants_psi_svd_beta ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_coefs" -> if (Ezfio.has_spindeterminants_psi_svd_coefs ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_alpha_numselected" -> if (Ezfio.has_spindeterminants_psi_svd_alpha_numselected ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_beta_numselected" -> if (Ezfio.has_spindeterminants_psi_svd_beta_numselected ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_alpha_numtoselect" -> if (Ezfio.has_spindeterminants_psi_svd_alpha_numtoselect ()) then "T" else "F"
 | "has_spindeterminants_psi_svd_beta_numtoselect" -> if (Ezfio.has_spindeterminants_psi_svd_beta_numtoselect ()) then "T" else "F"
 | "has_simulation_do_run" -> if (Ezfio.has_simulation_do_run ()) then "T" else "F"
 | "has_simulation_stop_time" -> if (Ezfio.has_simulation_stop_time ()) then "T" else "F"
 | "has_simulation_equilibration" -> if (Ezfio.has_simulation_equilibration ()) then "T" else "F"
 | "has_simulation_http_server" -> if (Ezfio.has_simulation_http_server ()) then "T" else "F"
 | "has_simulation_do_jast" -> if (Ezfio.has_simulation_do_jast ()) then "T" else "F"
 | "has_simulation_nucl_fitcusp_factor" -> if (Ezfio.has_simulation_nucl_fitcusp_factor ()) then "T" else "F"
 | "has_simulation_method" -> if (Ezfio.has_simulation_method ()) then "T" else "F"
 | "has_simulation_block_time" -> if (Ezfio.has_simulation_block_time ()) then "T" else "F"
 | "has_simulation_sampling" -> if (Ezfio.has_simulation_sampling ()) then "T" else "F"
 | "has_simulation_save_data" -> if (Ezfio.has_simulation_save_data ()) then "T" else "F"
 | "has_simulation_time_step" -> if (Ezfio.has_simulation_time_step ()) then "T" else "F"
 | "has_simulation_print_level" -> if (Ezfio.has_simulation_print_level ()) then "T" else "F"
 | "has_simulation_ci_threshold" -> if (Ezfio.has_simulation_ci_threshold ()) then "T" else "F"
 | "has_simulation_md5_key" -> if (Ezfio.has_simulation_md5_key ()) then "T" else "F"
 | "has_simulation_e_ref" -> if (Ezfio.has_simulation_e_ref ()) then "T" else "F"
 | "has_simulation_e_trial" -> if (Ezfio.has_simulation_e_trial ()) then "T" else "F"
 | "has_simulation_srmc_projection_time" -> if (Ezfio.has_simulation_srmc_projection_time ()) then "T" else "F"
 | "has_simulation_use_trexio" -> if (Ezfio.has_simulation_use_trexio ()) then "T" else "F"
 | "has_simulation_use_qmckl" -> if (Ezfio.has_simulation_use_qmckl ()) then "T" else "F"
 | "has_simulation_trexio_filename" -> if (Ezfio.has_simulation_trexio_filename ()) then "T" else "F"
 | "has_jastrow_jast_type" -> if (Ezfio.has_jastrow_jast_type ()) then "T" else "F"
 | "has_jastrow_jast_a_up_up" -> if (Ezfio.has_jastrow_jast_a_up_up ()) then "T" else "F"
 | "has_jastrow_jast_a_up_dn" -> if (Ezfio.has_jastrow_jast_a_up_dn ()) then "T" else "F"
 | "has_jastrow_jast_b_up_up" -> if (Ezfio.has_jastrow_jast_b_up_up ()) then "T" else "F"
 | "has_jastrow_jast_b_up_dn" -> if (Ezfio.has_jastrow_jast_b_up_dn ()) then "T" else "F"
 | "has_jastrow_mu_erf" -> if (Ezfio.has_jastrow_mu_erf ()) then "T" else "F"
 | "has_jastrow_jast_pen" -> if (Ezfio.has_jastrow_jast_pen ()) then "T" else "F"
 | "has_jastrow_jast_een_e_a" -> if (Ezfio.has_jastrow_jast_een_e_a ()) then "T" else "F"
 | "has_jastrow_jast_een_e_b" -> if (Ezfio.has_jastrow_jast_een_e_b ()) then "T" else "F"
 | "has_jastrow_jast_een_n" -> if (Ezfio.has_jastrow_jast_een_n ()) then "T" else "F"
 | "has_jastrow_jast_core_a1" -> if (Ezfio.has_jastrow_jast_core_a1 ()) then "T" else "F"
 | "has_jastrow_jast_core_a2" -> if (Ezfio.has_jastrow_jast_core_a2 ()) then "T" else "F"
 | "has_jastrow_jast_core_b1" -> if (Ezfio.has_jastrow_jast_core_b1 ()) then "T" else "F"
 | "has_jastrow_jast_core_b2" -> if (Ezfio.has_jastrow_jast_core_b2 ()) then "T" else "F"
 | "has_jastrow_jast_1b_type" -> if (Ezfio.has_jastrow_jast_1b_type ()) then "T" else "F"
 | "has_jastrow_jast_1btanh_pen" -> if (Ezfio.has_jastrow_jast_1btanh_pen ()) then "T" else "F"
 | "has_jastrow_jast_1berf_pen" -> if (Ezfio.has_jastrow_jast_1berf_pen ()) then "T" else "F"
 | "has_jastrow_jast_1bgauss_pen" -> if (Ezfio.has_jastrow_jast_1bgauss_pen ()) then "T" else "F"
 | "has_blocks_empty" -> if (Ezfio.has_blocks_empty ()) then "T" else "F"
 | "has_pseudo_ao_pseudo_grid" -> if (Ezfio.has_pseudo_ao_pseudo_grid ()) then "T" else "F"
 | "has_pseudo_do_pseudo" -> if (Ezfio.has_pseudo_do_pseudo ()) then "T" else "F"
 | "has_pseudo_mo_pseudo_grid" -> if (Ezfio.has_pseudo_mo_pseudo_grid ()) then "T" else "F"
 | "has_pseudo_pseudo_dz_k" -> if (Ezfio.has_pseudo_pseudo_dz_k ()) then "T" else "F"
 | "has_pseudo_pseudo_dz_kl" -> if (Ezfio.has_pseudo_pseudo_dz_kl ()) then "T" else "F"
 | "has_pseudo_pseudo_grid_rmax" -> if (Ezfio.has_pseudo_pseudo_grid_rmax ()) then "T" else "F"
 | "has_pseudo_pseudo_grid_size" -> if (Ezfio.has_pseudo_pseudo_grid_size ()) then "T" else "F"
 | "has_pseudo_pseudo_klocmax" -> if (Ezfio.has_pseudo_pseudo_klocmax ()) then "T" else "F"
 | "has_pseudo_pseudo_kmax" -> if (Ezfio.has_pseudo_pseudo_kmax ()) then "T" else "F"
 | "has_pseudo_pseudo_lmax" -> if (Ezfio.has_pseudo_pseudo_lmax ()) then "T" else "F"
 | "has_pseudo_pseudo_n_k" -> if (Ezfio.has_pseudo_pseudo_n_k ()) then "T" else "F"
 | "has_pseudo_pseudo_n_kl" -> if (Ezfio.has_pseudo_pseudo_n_kl ()) then "T" else "F"
 | "has_pseudo_pseudo_v_k" -> if (Ezfio.has_pseudo_pseudo_v_k ()) then "T" else "F"
 | "has_pseudo_pseudo_v_kl" -> if (Ezfio.has_pseudo_pseudo_v_kl ()) then "T" else "F"
 | x -> failwith (x^" : Unknown EZFIO function")
;;

let all_ezfio_messages = [ 
   "pseudo_pseudo_v_kl" ; 
   "pseudo_pseudo_v_k" ; 
   "pseudo_pseudo_n_kl" ; 
   "pseudo_pseudo_n_k" ; 
   "pseudo_pseudo_lmax" ; 
   "pseudo_pseudo_kmax" ; 
   "pseudo_pseudo_klocmax" ; 
   "pseudo_pseudo_grid_size" ; 
   "pseudo_pseudo_grid_rmax" ; 
   "pseudo_pseudo_dz_kl" ; 
   "pseudo_pseudo_dz_k" ; 
   "pseudo_mo_pseudo_grid" ; 
   "pseudo_do_pseudo" ; 
   "pseudo_ao_pseudo_grid" ; 
   "blocks_empty" ; 
   "jastrow_jast_1bgauss_pen" ; 
   "jastrow_jast_1berf_pen" ; 
   "jastrow_jast_1btanh_pen" ; 
   "jastrow_jast_1b_type" ; 
   "jastrow_jast_core_b2" ; 
   "jastrow_jast_core_b1" ; 
   "jastrow_jast_core_a2" ; 
   "jastrow_jast_core_a1" ; 
   "jastrow_jast_een_n" ; 
   "jastrow_jast_een_e_b" ; 
   "jastrow_jast_een_e_a" ; 
   "jastrow_jast_pen" ; 
   "jastrow_mu_erf" ; 
   "jastrow_jast_b_up_dn" ; 
   "jastrow_jast_b_up_up" ; 
   "jastrow_jast_a_up_dn" ; 
   "jastrow_jast_a_up_up" ; 
   "jastrow_jast_type" ; 
   "simulation_trexio_filename" ; 
   "simulation_use_qmckl" ; 
   "simulation_use_trexio" ; 
   "simulation_srmc_projection_time" ; 
   "simulation_e_trial" ; 
   "simulation_e_ref" ; 
   "simulation_md5_key" ; 
   "simulation_ci_threshold" ; 
   "simulation_print_level" ; 
   "simulation_time_step" ; 
   "simulation_save_data" ; 
   "simulation_sampling" ; 
   "simulation_block_time" ; 
   "simulation_method" ; 
   "simulation_nucl_fitcusp_factor" ; 
   "simulation_do_jast" ; 
   "simulation_http_server" ; 
   "simulation_equilibration" ; 
   "simulation_stop_time" ; 
   "simulation_do_run" ; 
   "spindeterminants_psi_svd_beta_numtoselect" ; 
   "spindeterminants_psi_svd_alpha_numtoselect" ; 
   "spindeterminants_psi_svd_beta_numselected" ; 
   "spindeterminants_psi_svd_alpha_numselected" ; 
   "spindeterminants_psi_svd_coefs" ; 
   "spindeterminants_psi_svd_beta" ; 
   "spindeterminants_psi_svd_alpha" ; 
   "spindeterminants_psi_svd_coefs_unique" ; 
   "spindeterminants_psi_svd_beta_unique" ; 
   "spindeterminants_psi_svd_alpha_unique" ; 
   "spindeterminants_n_svd_toselect" ; 
   "spindeterminants_n_svd_selected" ; 
   "spindeterminants_n_svd_coefs" ; 
   "spindeterminants_n_svd_coefs_unique" ; 
   "spindeterminants_psi_coef_matrix_values" ; 
   "spindeterminants_psi_coef_matrix_columns" ; 
   "spindeterminants_psi_coef_matrix_rows" ; 
   "spindeterminants_psi_det_beta" ; 
   "spindeterminants_psi_det_alpha" ; 
   "spindeterminants_n_states" ; 
   "spindeterminants_bit_kind" ; 
   "spindeterminants_n_int" ; 
   "spindeterminants_n_det" ; 
   "spindeterminants_n_det_beta" ; 
   "spindeterminants_n_det_alpha" ; 
   "nuclei_nucl_coord" ; 
   "nuclei_nucl_charge" ; 
   "nuclei_nucl_label" ; 
   "nuclei_nucl_num" ; 
   "electrons_elec_fitcusp_radius" ; 
   "electrons_elec_coord_pool_size" ; 
   "electrons_elec_coord_pool" ; 
   "electrons_elec_walk_num" ; 
   "electrons_elec_walk_num_tot" ; 
   "electrons_elec_beta_num" ; 
   "electrons_elec_alpha_num" ; 
   "mo_basis_mo_symmetry" ; 
   "mo_basis_mo_occ" ; 
   "mo_basis_mo_energy" ; 
   "mo_basis_mo_classif" ; 
   "mo_basis_mo_coef" ; 
   "mo_basis_mo_num" ; 
   "ao_basis_ao_expo" ; 
   "ao_basis_ao_coef" ; 
   "ao_basis_ao_power" ; 
   "ao_basis_ao_nucl" ; 
   "ao_basis_ao_prim_num" ; 
   "ao_basis_ao_num" ; 
   "properties_ci_dress_mu_opt" ; 
   "properties_ci_h_matrix_diag" ; 
   "properties_ci_h_matrix" ; 
   "properties_ci_overlap_matrix" ; 
   "properties_ci_h_psidet" ; 
   "properties_ci_overlap_psidet" ; 
   "properties_ci_dress_mu" ; 
   "properties_ci_dress" ; 
   "properties_psi_norm" ; 
   "properties_emudiff" ; 
   "properties_ci_h_transcor_psi" ; 
   "properties_drift_mod" ; 
   "properties_pop_weight" ; 
   "properties_wf_extension" ; 
   "properties_dipole" ; 
   "properties_e_loc_zv" ; 
   "properties_e_loc" ; 
   "properties_e_kin" ; 
   "properties_e_pot" ; 
   "properties_e_nucl" ; 
   "ezfio_last_library" ; 
   "ezfio_library" ; 
   "ezfio_user" ; 
   "ezfio_creation" ; 
]
