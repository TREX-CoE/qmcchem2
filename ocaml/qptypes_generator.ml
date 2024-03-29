let global_replace x =
  x
  |> Str.global_replace (Str.regexp "Float.to_string") "string_of_float"
  |> Str.global_replace (Str.regexp "Float.of_string") "float_of_string"
  |> Str.global_replace (Str.regexp "Int.to_string") "string_of_int"
  |> Str.global_replace (Str.regexp "Int.of_string") "int_of_string"
  |> Str.global_replace (Str.regexp "Int.to_bytes") "bytes_of_int"
  |> Str.global_replace (Str.regexp "Int64.to_bytes") "bytes_of_int64"
  |> Str.global_replace (Str.regexp "Float.to_bytes") "bytes_of_float"
  |> Str.global_replace (Str.regexp "Float.of_bytes") "float_of_bytes"
  |> Str.global_replace (Str.regexp "Int.of_bytes") "int_of_bytes"
  |> Str.global_replace (Str.regexp "Int64.of_bytes") "int64_of_bytes"
  |> Str.global_replace (Str.regexp "String.\\(to\\|of\\)_string") ""
  |> Str.global_replace (Str.regexp "String.to_bytes") "Bytes.of_string"
  |> Str.global_replace (Str.regexp "String.of_bytes") "Bytes.to_string"

let input_data = "
* Positive_float : float
  if not (x >= 0.) then
    raise (Invalid_argument (Printf.sprintf \"Positive_float : (x >= 0.) : x=%f\"  x));

* Strictly_positive_float : float
  if not (x > 0.) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_positive_float : (x > 0.) : x=%f\" x));

* Negative_float : float
  if not (x <= 0.) then
    raise (Invalid_argument (Printf.sprintf \"Negative_float : (x <= 0.) : x=%f\" x));

* Strictly_negative_float : float
  if not (x < 0.) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_negative_float : (x < 0.) : x=%f\" x));

* Positive_int64 : int64
  if not (x >= 0L) then
    raise (Invalid_argument (Printf.sprintf \"Positive_int64 : (x >= 0L) : x=%s\" (Int64.to_string x)));

* Positive_int : int
  if not (x >= 0) then
    raise (Invalid_argument (Printf.sprintf \"Positive_int : (x >= 0) : x=%d\" x));

* Strictly_positive_int : int
  if not (x > 0) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_positive_int : (x > 0) : x=%d\" x));

* Negative_int : int
  if not (x <= 0) then
    raise (Invalid_argument (Printf.sprintf \"Negative_int : (x <= 0) : x=%d\" x));
  assert (x <= 0) ;

* Det_coef : float
  if (x < -1.) || (x > 1.) then
    raise (Invalid_argument (Printf.sprintf \"Det_coef : (-1. <= x <= 1.) : x=%f\" x));

* Normalized_float : float
  if (x < 0.) || (x > 1.) then
    raise (Invalid_argument (Printf.sprintf \"Normalized_float : (0. <= x <= 1.) : x=%f\" x));

* Strictly_negative_int : int
  if not (x < 0) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_negative_int : (x < 0) : x=%d\" x));

* Non_empty_string : string
  if (x = \"\") then
    raise (Invalid_argument \"Non_empty_string\");


* Det_number_max : int
  assert (x > 0) ;
  if (x > 100_000_000) then
    warning \"More than 100 million determinants\";

* States_number : int
  assert (x > 0) ;
  if (x > 1000) then
    warning \"More than 1000 states\";

* Bit_kind_size : int
  begin match x with
  | 8 | 16 | 32 | 64 -> ()
  | _ -> raise (Invalid_argument \"Bit_kind_size should be (8|16|32|64).\")
  end;

* Bit_kind : int
  begin match x with
  | 1 | 2 | 4 | 8 -> ()
  | _ -> raise (Invalid_argument \"Bit_kind should be (1|2|4|8).\")
  end;

* Bitmask_number : int
  assert (x > 0) ;

* MO_coef : float

* MO_occ : float
  if x < 0. then 0.  else
  if x > 2. then 2.  else

* AO_coef : float

* AO_expo : float
  if (x < 0.) then
    raise (Invalid_argument (Printf.sprintf \"AO_expo : (x >= 0.) : x=%f\" x));

* AO_prim_number : int
  assert (x > 0) ;

* Threshold : float
  assert (x >= 0.) ;
  assert (x <= 1.) ;

* PT2_energy : float
  assert (x >=0.) ;

* Elec_alpha_number : int
  assert (x > 0) ;

* Elec_beta_number : int
  assert (x >= 0) ;

* Elec_number : int
  assert (x > 0) ;

* MD5 : string
  assert ((String.length x) = 32);
  assert (
    let a =
      Array.init (String.length x) (fun i -> x.[i])
    in
    Array.fold_left (fun accu x -> accu && (x < 'g')) true a
    );

* Rst_string : string


* Weight : float
  assert (x >= 0.) ;

* Block_id : int
  assert (x > 0) ;

* Compute_node : string
  assert (x <> \"\") ;

"


let input_ezfio = "
* MO_number : int
  mo_basis_mo_num
  1 : 10_000
  More than 10_000 MOs

* AO_number : int
  ao_basis_ao_num
  1 : 10_000
  More than 10_000 AOs

* Nucl_number : int
  nuclei_nucl_num
  1 : 10_000
  More than 10_000 nuclei

* N_int_number : int
  spindeterminants_n_int
  1 : 30
  N_int > 30

* Det_number : int
  spindeterminants_n_det
  1 : 100_000_000
  More than 100 million determinants

"


let untouched = "
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

"

let template = format_of_string "
module %s : sig
  type t [@@deriving sexp]
  val to_%s : t -> %s
  val of_%s : %s %s -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
  val of_bytes  : bytes -> t
end = struct
  type t = %s [@@deriving sexp]
  let to_%s x = x
  let of_%s %s x = ( %s x )
  let to_string x = %s.to_string x
  let to_bytes  x = %s.to_bytes  x
  let of_bytes  b = %s.of_bytes  b
end

"



let parse_input input=
  print_string "open Sexplib.Std\nlet warning = print_string\n" ;
  let rec parse result = function
    | [] -> result
    | ( "" , ""   )::tail -> parse result tail
    | ( t  , text )::tail ->
        let name,typ,params,params_val =
          match String.split_on_char ':' t with
          | [name;typ] -> (name,typ,"","")
          | name::typ::params::params_val -> (name,typ,params,
            (String.concat ":" params_val) )
          | _ -> assert false
        in
        let typ  = String_ext.strip typ
        and name = String_ext.strip name in
        let typ_cap = String.capitalize_ascii typ in
        let newstring = Printf.sprintf template name typ typ typ params_val typ typ
          typ typ params ( String_ext.strip text ) typ_cap typ_cap typ_cap
        in
        List.rev (parse (newstring::result) tail )
  in
     String_ext.split ~on:'*' input
  |> List.map (String_ext.lsplit2_exn ~on:'\n')
  |> parse []
  |> String.concat  ""
  |> global_replace
  |> print_string



let ezfio_template = format_of_string "
module %s : sig
  type t [@@deriving sexp]
  val to_%s : t -> %s
  val get_max : unit -> %s
  val of_%s : ?min:%s -> ?max:%s -> %s -> t
  val to_string : t -> string
  val to_bytes  : t -> bytes
end = struct
  type t = %s [@@deriving sexp]
  let to_string x = %s.to_string x
  let to_bytes  x = %s.to_bytes  x
  let get_max () =
    if (Ezfio.has_%s ()) then
      Ezfio.get_%s ()
    else
      %s
  let get_min () =
      %s
  let to_%s x = x
  let of_%s ?(min=get_min ()) ?(max=get_max ()) x =
    begin
      assert (x >= min) ;
      if (x > %s) then
        warning \"%s\";
      begin
        match max with
        | %s -> ()
        | i  ->
          if ( x > i ) then
            raise (Invalid_argument (Printf.sprintf \"%s: %%s\" (%s.to_string x) ))
      end ;
      x
    end
end
"

(*
  val of_bytes  : bytes -> t
  let of_bytes  x = %s.of_bytes  x
*)

let parse_input_ezfio input=
  let parse s =
    match (
      String_ext.split s ~on:'\n'
      |> List.filter (fun x -> (String_ext.strip x) <> "")
    ) with
    | [] -> ""
    | a :: b :: c :: d :: [] ->
      begin
        let (name,typ) = String_ext.lsplit2_exn ~on:':' a
        and ezfio_func = b
        and (min, max) = String_ext.lsplit2_exn ~on:':' c
        and msg = d
        in
        let (name, typ, ezfio_func, min, max, msg) =
        match List.map String_ext.strip [ name ; typ ; ezfio_func ; min ; max ; msg ] with
        | [ name ; typ ; ezfio_func ; min ; max ; msg ] -> (name, typ, ezfio_func, min, max, msg)
        | _ -> assert false
        in
        let typ_cap = String.capitalize_ascii typ in
        Printf.sprintf ezfio_template
          name typ typ typ typ typ typ typ typ typ_cap typ_cap
          ezfio_func ezfio_func max min typ typ max msg min name typ_cap
      end
    | _ -> failwith "Error in input_ezfio"
  in
     String_ext.split ~on:'*' input
  |> List.map parse
  |> String.concat ""
  |> global_replace
  |> print_string



(** EZFIO *)
let input_lines filename =
  let ic = open_in filename in
  let result = String_ext.input_lines ic in
  close_in ic;
  result


let create_ezfio_handler () =
  let lines =
    input_lines "ezfio.ml"
    (* /!\ Change when ezfio.ml changes *)
    |> List.mapi (fun i l -> if i > 441 then Some l else None)
    |> List.filter (fun x -> x <> None)
    |> List.map (fun x ->
        match x with
        | Some x -> x
        | None -> assert false)
  in
  let functions =
    List.map (fun x ->
      match String.split_on_char ' ' x with
      | _ :: x :: "()" :: "=" :: f :: dir :: item :: _-> (x, f, dir, item)
      | _ :: x ::         "=" :: f :: dir :: item :: _-> (x, f, dir, item)
      | _ -> ("","","","")
    ) lines
  in
  let has_functions =
    List.filter (fun (x,_,_,_) -> String.sub x 0 4 = "has_") functions
  and get_functions =
    List.filter (fun (x,_,_,_) -> String.sub x 0 4 = "get_") functions
  in
  let chop s =
    match (Str.split_delim (Str.regexp ";;") s) with
    | x :: _ -> x
    | _ -> assert false
  in

  let result =
    [ "let decode_ezfio_message msg =
match msg with " ] @
    (
      List.map (fun (x,f,d,i) ->
        let i = chop i in
        if (String.sub f ((String.length f)-6) 6 = "_array") then
          Printf.sprintf " | \"%s\" ->
             Ezfio.read_string_array %s %s
             |> Ezfio.flattened_ezfio
             |> Array.to_list
             |> String.concat \" \"" x d i
        else
          Printf.sprintf " | \"%s\" -> Ezfio.read_string %s %s" x d i
      ) get_functions
    ) @ (
      List.map (fun (x,_,_,_) ->
        Printf.sprintf " | \"%s\" -> if (Ezfio.%s ()) then \"T\" else \"F\"" x x
      ) has_functions
    )
    @ [" | x -> failwith (x^\" : Unknown EZFIO function\")\n;;"  ;
       "" ; "let all_ezfio_messages = [ " ] @
    (
      List.rev_map (fun (x,_,_,_) ->
        Printf.sprintf "   \"%s\" ; " (String.sub x 4 ((String.length x)-4))
      ) has_functions
    ) @ ["]"]
  in
  String.concat "\n" result
  |> print_endline

(** Main *)
let () =
  print_endline untouched;
  parse_input input_data ;
  parse_input_ezfio input_ezfio;
  create_ezfio_handler ()




