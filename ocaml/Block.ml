open Qptypes

type t =
  { property       : Property.t ;
    value          : Sample.t        ;
    weight         : Weight.t        ;
    compute_node   : Compute_node.t  ;
    pid            : int             ;
    block_id       : Block_id.t      ;
  }

let re =
  Str.regexp "[ |#|\n]+"

let of_string s =

   try
     let lst =
        Str.split re s
        |> List.rev
     in
      match lst with
      | b :: pid :: c:: p :: w :: v :: [] -> Some
        { property = Property.of_string p ;
          value    = Sample.of_float (float_of_string v) ;
          weight   = Weight.of_float (float_of_string w) ;
          compute_node = Compute_node.of_string c;
          pid          = int_of_string pid;
          block_id     = Block_id.of_int (int_of_string b) ;
        }
      | b :: pid :: c:: p :: w :: v ->
          let v =
            List.rev v
            |> Array.of_list
            |> Array.map float_of_string
          in
          let dim =
            Array.length v
          in
        Some
        { property = Property.of_string p ;
          value    = Sample.of_float_array ~dim v ;
          weight   = Weight.of_float (float_of_string w) ;
          compute_node = Compute_node.of_string c;
          pid          = int_of_string pid;
          block_id     = Block_id.of_int (int_of_string b) ;
        }
      | _ -> None
   with
   | _ -> None



let to_short_string b =
    Printf.sprintf "%s # %s %d %d"
      (Property.to_string     b.property)
      (Compute_node.to_string b.compute_node)
      b.pid
      (Block_id.to_int        b.block_id)


let to_string b =
    Printf.sprintf "%s %s # %s %s %s %d"
      (Sample.to_string       b.value )
      (Weight.to_float        b.weight       |> string_of_float)
      (Property.to_string     b.property)
      (Compute_node.to_string b.compute_node)
      (string_of_int          b.pid)
      (Block_id.to_int        b.block_id)




let zero =
  bytes_of_int 0

let to_bytes b =
  (* [ Length of b
       [ Length of value ;
         Value ;
         Length of weight ;
         Weight ;
         ...  ] ]  *)
  let l =
    [  Property.to_bytes      b.property ;
       Sample.to_bytes        b.value  ;
       Weight.to_bytes        b.weight ;
       bytes_of_int           b.pid ;
       Block_id.to_bytes      b.block_id ;
       Compute_node.to_bytes  b.compute_node ]
    |> List.map (fun x -> [ bytes_of_int (Bytes.length x) ;  x ] )
    |> List.concat
  in
  let result =
    Bytes.concat Bytes.empty (zero :: l)
  in
  Bytes.set_int64_ne result 0 (Int64.of_int ((Bytes.length result) - 8));
  result


let read_bytes b idx =
  (* Reads m, the first 8 bytes as an int64 containing the number of bytes to read.
     Then, read the next m bytes and return a tuple containing the decoded data and the rest.
  *)
  let l = (Bytes.length b) - idx in
  if l < 8 then
     None
  else
     let m =
       Bytes.get_int64_ne b idx
       |> Int64.to_int
     in
     try
       Some (Bytes.sub b (idx+8) m, idx+8+m)
     with Invalid_argument _ -> None



let of_bytes ?(idx=0) b =
  let get_x s idx =
    match read_bytes s idx with
    | Some ( data, i1) -> data, i1
    | _ -> raise Exit
  in


  let result =
    try
      let property, idx = get_x b idx in
      let value   , idx = get_x b idx in
      let weight  , idx = get_x b idx in
      let pid     , idx = get_x b idx in
      let block_id, idx = get_x b idx in
      let compute_node, i5 = get_x b idx in
      Some
          { property = Property.of_bytes property;
            value    = Sample.of_bytes value;
            weight   = Weight.of_bytes weight;
            pid          = int_of_bytes pid;
            block_id     = Block_id.of_bytes block_id;
            compute_node = Compute_node.of_bytes compute_node;
          }
    with Exit -> None
  in
  result


let of_string_or_bytes s =
  if Qmcchem_config.binary_io then
     Bytes.of_string s
     |> of_bytes
  else
     of_string s

let dir_name = lazy(
  let ezfio_filename =
     Lazy.force Qputils.ezfio_filename
  in
  let md5 =
     QmcMd5.hash ()
  in
  let d = Filename.concat ezfio_filename "blocks" in
  if not ( Sys.file_exists d ) then
    Unix.mkdir d 0o755;
  List.fold_right Filename.concat
    [ ezfio_filename ; "blocks" ; md5 ; Filename.dir_sep ] ""
)


(* Fetch raw data from the EZFIO file *)
let _raw_data =
  ref None


let update_raw_data ?(locked=true) () =
  (* Create array of files to read *)
  let dir_name =
    Lazy.force dir_name
  in
  let files =
     let result =
       if Sys.file_exists dir_name && Sys.is_directory dir_name then
         begin
           Sys.readdir dir_name
           |> Array.map (fun x -> dir_name^x)
           |> Array.to_list
         end
       else []
     in
     if locked then
       result
     else
       List.filter (fun x ->
         try
          let _ =
            Str.search_backward (Str.regexp "locked") x ((String.length x) - 1)
          in false
         with
         | Not_found -> true
       ) result
  in

  if Qmcchem_config.binary_io then
    begin
      let result =
        let rec aux buf idx accu =
          (* Read one block *)
          match read_bytes buf idx with
          | None -> List.rev accu
          | Some (item, new_idx) ->
              match of_bytes item with
              | None   -> List.rev accu
              | Some item -> (aux [@tailcall]) buf new_idx (item::accu)
        in
        List.concat_map (fun filename ->
          let ic = open_in filename in
          let length = in_channel_length ic in
          let result =
            if length > 0 then
                let buf = Bytes.create length in
                really_input ic buf 0 length;
                aux buf 0 []
            else []
          in
          close_in ic;
          result ) files
      in
      result
    end
  else
    begin
      let rec transform new_list = function
      | [] -> new_list
      | head :: tail ->
        let head = String.trim head in
        let item = of_string head in
        match item with
          | None   -> transform new_list tail
          | Some x -> transform (x::new_list) tail
      in

      let result =
        let rec aux ic accu =
          let l =
            try
              Some (input_line ic)
            with
            | End_of_file -> None
          in
          match l with
          | None -> List.rev accu
          | Some l -> (aux [@tailcall]) ic (l::accu)
          in
          List.concat_map (fun filename ->
            let ic = open_in filename in
            let result = aux ic [] in
            close_in ic;
            result ) files
        |> transform []
      in
      result
    end


let to_string_or_bytes b =
  if Qmcchem_config.binary_io then
    to_bytes b
    |> Bytes.to_string
  else
    to_string b


let raw_data ?(locked=true) () =
  match !_raw_data with
  | Some x -> x
  | None   ->
    let result =
      update_raw_data ~locked ()
    in
    _raw_data := Some result;
    result



let properties = lazy (
  let h = Hashtbl.create 63 in
  List.iter (fun x ->
    Hashtbl.replace h (Property.to_string x.property) x.property)
    (raw_data ());
  Hashtbl.fold (fun k v a -> v :: a) h []
)


