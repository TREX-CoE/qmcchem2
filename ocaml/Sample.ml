open Sexplib.Std

type t =
  | One_dimensional  of float
  | Multidimensional of (float array * int)
[@@deriving sexp]

let dimension = function
  | One_dimensional _ -> 1
  | Multidimensional (_,d) -> d

let to_float ?idx x =
  match (idx,x) with
  | None  , One_dimensional x
  | Some 0, One_dimensional x  -> x
  | Some i, One_dimensional x  ->
      failwith "Index should not be specified in One_dimensional"
  | None  , Multidimensional (x,_) -> x.(0)
  | Some i, Multidimensional (x,s) when i < s  -> x.(i)
  | Some i, Multidimensional (x,s) ->
      Printf.sprintf "Index out of bounds in Multidimensional
      %d not in [0,%d[ " i s
      |> failwith

let to_float_array = function
  | One_dimensional _ -> failwith "Should be Multidimensional"
  | Multidimensional (x,_) -> x

let of_float x =
  One_dimensional x

let of_float_array ~dim x =
  if (Array.length x) <> dim then
    failwith "Inconsistent array size in of_float_array"
  else
    match dim with
    | 1 -> One_dimensional x.(0)
    | _ -> Multidimensional (x, dim)

let to_string = function
  | One_dimensional  x -> string_of_float x
  | Multidimensional (x,_) ->
      Array.map string_of_float x
      |> Array.to_list
      |> String.concat " "

let to_bytes = function
  | One_dimensional  x -> Qptypes.bytes_of_float x
  | Multidimensional (x,_) ->
      let b = Bytes.create (8 * Array.length x) in
      Array.iteri (fun i x ->
        Int64.bits_of_float x
        |> Bytes.set_int64_ne b (i*8) ) x;
      b

let of_bytes b =
  match Bytes.length b with
  | 8 -> let x = Qptypes.float_of_bytes b in
         One_dimensional x
  | l -> let len = l/8 in
let result = 
         Multidimensional ( Array.init len (fun i ->
                               Bytes.get_int64_ne b (i*8)
                               |> Int64.float_of_bits ),
                           len )
in
result
