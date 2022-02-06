open Qptypes

type t =
  | Property    of Block.t
  | Walkers     of Compute_node.t * int * (float array) array
  | Register    of Compute_node.t * int
  | Unregister  of Compute_node.t * int
  | Test
  | GetWalkers  of Strictly_positive_int.t
  | Ezfio       of string
  | Error       of string


let of_string_list m =
  try
    match m with
    | [ "cpu"        ; c ; pid ; b ; "1" ; v ] ->
        let open Block in
            Property
            { property     = Property.Cpu;
              value        = Sample.of_float (float_of_string v) ;
              weight       = Weight.of_float 1.;
              compute_node = Compute_node.of_string c;
              pid          = int_of_string pid;
              block_id     = Block_id.of_int (int_of_string b);
            }
    | [ "accep" ; c ; pid ; b ; "1" ; v ] ->
        let open Block in
            Property
            { property     = Property.Accep;
              value        = Sample.of_float (float_of_string v) ;
              weight       = Weight.of_float 1.;
              compute_node = Compute_node.of_string c;
              pid          = int_of_string pid;
              block_id     = Block_id.of_int (int_of_string b);
            }
    | [ prop ; c ; pid ; b ; w ; v ] ->
        let open Block in
            Property
            { property     = Property.of_string prop;
              value        = Sample.of_float (float_of_string v);
              weight       = Weight.of_float (float_of_string w);
              compute_node = Compute_node.of_string c;
              pid          = int_of_string pid;
              block_id     = Block_id.of_int (int_of_string b);
            }
    | "elec_coord" :: c :: pid :: _ :: n ::walkers ->
      begin
        let elec_num =
          Lazy.force Qputils.elec_num
        and n =
          int_of_string n
        in
        assert (n = List.length walkers);
        let rec build_walker accu = function
        | (0,tail) ->
            let result =
              List.rev accu
              |> List.rev_map float_of_string
              |> List.rev
              |> Array.of_list
            in
            (result, tail)
        | (n,head::tail)  ->
            build_walker (head::accu) (n-1, tail)
        | _ -> failwith "Bad walkers"
        in
        let rec build accu = function
        | [] -> Array.of_list accu
        | w ->
          let (result, tail) =
            build_walker [] (3*elec_num+3, w)
          in
          build (result::accu) tail
        in
        Walkers (Compute_node.of_string c, int_of_string pid, build [] walkers)
      end
    | [ "get_walkers" ; n ] -> GetWalkers (n |> int_of_string |> Strictly_positive_int.of_int)
    | [ "register"   ; c ; pid  ] -> Register   (Compute_node.of_string c, int_of_string pid)
    | [ "unregister" ; c ; pid  ] -> Unregister (Compute_node.of_string c, int_of_string pid)
    | [ "Test" ] -> Test
    | [ "Ezfio" ; ezfio_msg ] -> Ezfio ezfio_msg
    | prop :: c :: pid :: b :: d :: w :: "bin" :: block :: [] ->
        (* Block in binary format *)
        let property =
          Property.of_string prop
        in
        begin
          assert (not (Property.is_scalar property));
          match Block.of_bytes ~idx:8 (Bytes.unsafe_of_string block)  with
          | Some block -> Property block
          | None -> failwith "Invalid block"
        end
    | prop :: c :: pid :: b :: d :: w :: l ->
        (* Bock in text format *)
        let property =
          Property.of_string prop
        in
        begin
          assert (not (Property.is_scalar property));
          let a =
            Array.of_list l
            |> Array.map float_of_string
          and dim =
            int_of_string d
          in
          assert (Array.length a = dim);
          let open Block in
            Property
            { property     = property ;
              value        = Sample.of_float_array ~dim a;
              weight       = Weight.of_float (float_of_string w);
              compute_node = Compute_node.of_string c;
              pid          = int_of_string pid;
              block_id     = Block_id.of_int (int_of_string b);
            }
        end
    | l -> Error (String.concat ":" l)
  with
  | Assert_failure (l,_,_) -> Error l
  | _ -> Error "Unknown error"



let to_string = function
  | Property   b         -> "Property   : "^(Block.to_string b)
  | Walkers    (h,p,w)   -> Printf.sprintf "Walkers      : %s %d : %d walkers"
                              (Compute_node.to_string h) p (Array.length w)
  | GetWalkers n         -> Printf.sprintf "GetWalkers %d" (Strictly_positive_int.to_int n)
  | Register   (h,p)     -> Printf.sprintf "Register   : %s %d"
                              (Compute_node.to_string h) p
  | Unregister (h,p)     -> Printf.sprintf "Unregister   : %s %d"
                              (Compute_node.to_string h) p
  | Test                 -> "Test"
  | Ezfio msg            -> "Ezfio "^msg
  | Error msg            -> "Error "^msg


let create m =
   of_string_list m

