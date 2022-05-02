open Qptypes

type t =
    { property  :  Property.t ;
      data      :  Block.t list;
    }


module Average = struct
  include Sample
end

module Error = struct
  include Sample
end

module Variance = struct
  include Sample
end

module Skewness: sig
  type t
  val to_float : t -> float
  val of_float : float -> t
  val to_string : t -> string
end = struct
  type t = float
  let to_string = string_of_float
  let to_float x = x
  let of_float x = x
end

module Kurtosis: sig
  type t
  val to_float : t -> float
  val of_float : float -> t
  val to_string : t -> string
end = struct
  type t = float
  let to_string = string_of_float
  let to_float x = x
  let of_float x = x
end

module GaussianDist: sig
  type t
  val create : mu:Average.t -> sigma2:Variance.t -> t
  val eval : g:t -> x:float -> float
end = struct
  type t = { mu: Average.t ; sigma2: Variance.t }
  let create ~mu ~sigma2 =
    { mu ; sigma2 }
  let eval ~g ~x =
    let { mu ; sigma2 } =
      g
    in
    let mu =
      Average.to_float mu
    and sigma2 =
      Variance.to_float sigma2
    in
    let x2 =
      (x -. mu) *. ( x -. mu) /. sigma2
    in
    let pi =
      acos (-1.)
    in
    let c =
      1. /. (sqrt (sigma2 *. (pi +. pi)))
    in
    c *. exp ( -0.5 *. x2)

end


let hashtbl_to_alist table =
  Hashtbl.fold (fun k v a ->  (k,v) :: a) table []

let hashtbl_change table key f =
  let elt =
    try
      Some (Hashtbl.find table key)
    with
    | Not_found -> None
  in
  let new_elt = f elt in
  match new_elt with
    | None -> Hashtbl.remove table key
    | Some value -> Hashtbl.replace table key value


(** Compute average *)
let average { property ; data } =
  if Property.is_scalar property then
    let (num,denom) =
      List.fold_left (fun (an, ad) x ->
        let num =
          (Weight.to_float x.Block.weight) *. (Sample.to_float x.Block.value)
        and den =
          (Weight.to_float x.Block.weight)
        in (an +. num, ad +. den)
      ) (0., 0.) data
    in
    num /. denom
    |> Average.of_float
  else
    let dim =
      match data with
      | [] -> 1
      | x :: tl -> Sample.dimension x.Block.value
    in
    let (num,denom) =
      List.fold_left (fun (an, ad) x ->
        let num =
          Array.map (fun y -> (Weight.to_float x.Block.weight) *. y)
            (Sample.to_float_array x.Block.value)
        and den = (Weight.to_float x.Block.weight)
        in ( Array.mapi (fun i y -> y +. num.(i)) an , ad +. den)
      ) (Array.make dim 0. , 0.) data
    in
    let denom_inv =
      1. /. denom
    in
    Array.map (fun x -> x *. denom_inv) num
    |> Average.of_float_array ~dim



(** Compute variance *)
let variance { property ; data } =
  let average = average { property ; data } in
  if Property.is_scalar property then
    let average = Average.to_float average in
    let (num,denom) =
      List.fold_left (fun (an, ad) x ->
        let num =
          let z = (Sample.to_float x.Block.value) -. average in
          (Weight.to_float x.Block.weight) *. z *. z
        and den =
          (Weight.to_float x.Block.weight)
        in (an +. num, ad +. den)
      ) (0., 0.) data
    in
    num /. denom
    |> Variance.of_float
  else
    let average = Average.to_float_array average in
    let dim =
      match data with
      | [] -> 1
      | x :: tl -> Sample.dimension x.Block.value
    in
    let (num,denom) =
      List.fold_left (fun (an, ad) x ->
        let num =
          Array.mapi (fun i y ->
             let z = (y -. average.(i)) in
              z *. z
          )
          (Sample.to_float_array x.Block.value)
        and den = (Weight.to_float x.Block.weight)
        in ( Array.mapi (fun i y -> y +. num.(i)) an , ad +. den)
      ) (Array.make dim 0. , 0.) data
    in
    let denom_inv =
      1. /. denom
    in
    Array.map (fun x -> x *. denom_inv) num
    |> Variance.of_float_array ~dim





(** Compute sum (for CPU/Wall time) *)
let sum { property ; data } =
  List.fold_left (fun accu x ->
       let num = (Weight.to_float x.Block.weight) *. (Sample.to_float x.Block.value)
       in accu +. num
     ) 0. data

(** Compute sum of the weights *)
let total_weight { property ; data } =
  List.fold_left (fun accu x ->
       let num = (Weight.to_float x.Block.weight)
       in accu +. num
     ) 0. data



(** Calculation of the average and error bar *)
let ave_error { property ; data } =

  (* sum:   \sum_k  x_k *. w_k
     ansum: \sum_k  w_k
     avsum: \sum_k  x_k *. w_k
     avcu0: avsum / ansum
     avsq:  \sum_k (1. -. (w_k /. ansum_k)) *. (x_k -. avcu0)^2 *. w_k)
  *)
  let rec loop ~sum ~avsq ~ansum ~avsum ~n ?idx = function
    | [] ->
      begin
        if (n > 0.) then
           ( Average.of_float (sum /. ansum),
             Some (Error.of_float (sqrt ( abs_float ( avsq /.( ansum *. n)))) ))
        else
          ( Average.of_float (sum /. ansum), None)
      end
    | (x, w) :: tail ->
      begin
        let avcu0 = avsum /. ansum in
        let xw = x *. w in
        let ansum, avsum, sum =
          ansum +. w ,
          avsum +. xw ,
          sum   +. xw
        in
        loop tail
          ~sum:sum
          ~avsq:(avsq +. (1. -. (w /. ansum)) *. (x -. avcu0)
                  *. (x -. avcu0) *. w)
          ~avsum:avsum
          ~ansum:ansum
          ~n:(n +. 1.)
      end
  in

  let ave_error_scalar = function
    | [] -> (Average.of_float 0., None)
    | (x,w) :: tail ->
        loop tail
            ~sum:(x *. w)
            ~avsq:0.
            ~ansum:w
            ~avsum:(x *. w)
            ~n:0.
  in

  if (Property.is_scalar property) then
    List.rev_map (fun x ->
      (Sample.to_float  x.Block.value,
       Weight.to_float  x.Block.weight)
    ) data
    |> ave_error_scalar
  else
    match data with
    | [] -> (Average.of_float 0., None)
    | head::tail as list_of_samples ->
      let dim =
        head.Block.value
        |> Sample.dimension
      in
      let result =
        Array.init dim (fun idx ->
          List.rev_map (fun x ->
            (Sample.to_float  ~idx  x.Block.value,
             Weight.to_float        x.Block.weight)
          ) list_of_samples
          |> ave_error_scalar
        )
      in
      ( Array.map (fun (x,_) -> Average.to_float x) result
        |> Average.of_float_array ~dim ,
        if (Array.length result < 2) then
          None
        else
          Some (Array.map (function
            | (_,Some y) -> Error.to_float y
            | (_,None)   -> 0.) result
          |> Average.of_float_array ~dim)
      )

(** Computes the first 4 centered cumulants (zero mean) *)
let centered_cumulants { property ; data } =
  let ave =
    average { property ; data }
    |> Average.to_float
  in
  let centered_data =
     List.rev_map (fun x ->
      ( (Weight.to_float x.Block.weight),
        (Sample.to_float x.Block.value) -. ave )
      ) data
     |> List.rev
  in
  let var  =
     let (num, denom) =
     List.fold_left (fun (a2, ad) (w,x) ->
       let x2 = x *. x
       in
       let var = w *. x2
       and den = w
       in (a2 +. var, ad +. den)
     ) (0., 0.) centered_data
     in num /. denom
  in
  let centered_data =
     let sigma_inv =
       1. /. (sqrt var)
     in
     List.rev_map (fun x ->
      ( (Weight.to_float x.Block.weight),
        ( (Sample.to_float x.Block.value) -. ave ) *. sigma_inv )
      ) data
     |> List.rev
  in
  let (cum3,cum4) =
     let (cum3, cum4, denom) =
     List.fold_left (fun (a3, a4, ad) (w,x) ->
       let x2 = x *. x
       in
       let cum3 = w *. x2 *. x
       and cum4 = w *. x2 *. x2
       and den = w
       in (a3 +. cum3, a4 +. cum4, ad +. den)
     ) (0., 0., 0.) centered_data
     in
     ( cum3 /. denom, cum4 /. denom -. 3. )
  in
  [| ave ; var  ; cum3 ;  cum4 |]


(** Build from raw data. Range values are given in percent. *)
let of_raw_data ?(locked=true) ~range ~clean property =
    let raw_data =
      let data =
         Block.raw_data ~locked ()
         |> List.filter (fun x -> x.Block.property = property)
      in
      match clean with
      | None -> data
      | Some clean ->
          let average =
            average { property ; data }
          in
          let variance =
            variance { property ; data }
          in
          let result =
            if Property.is_scalar property then
              let var = Variance.to_float variance in
              let sigma = sqrt var in
              let ave = Average.to_float average in
                List.filter (fun x ->
                          abs_float ((Sample.to_float x.Block.value) -. ave) < clean *. sigma
                ) data
            else
              let sigma =
                  Array.map sqrt (Variance.to_float_array variance)
              in
              let ave = Average.to_float_array average in
              List.filter (fun x ->
                Array.mapi (fun i y ->
                  abs_float (y -. ave.(i)) < clean *. sigma.(i)
                ) (Sample.to_float_array x.Block.value)
                |> Array.fold_left (fun x y -> x && y) true
              ) data
          in
          match result with
          | [] -> data
          | _ -> result
    in
    let data =
         raw_data
         |> List.sort (fun x y ->
              let x = Block_id.to_int x.Block.block_id in
              let y = Block_id.to_int y.Block.block_id in
              if (x>y) then 1
              else if (x<y) then -1
              else 0)
    in

    let data_in_range rmin rmax =

        let total_weight =
          List.fold_left (fun accu x ->
              (Weight.to_float x.Block.weight) +. accu
            ) 0.  data
        in

        let wmin, wmax =
          rmin *. total_weight *. 0.01,
          rmax *. total_weight *. 0.01
        in

        let (_, new_data) =
            List.fold_left (fun (wsum, l) x ->
                if (wsum > wmax) then
                  (wsum,l)
                else
                  begin
                    let wsum_new =
                        wsum +. (Weight.to_float x.Block.weight)
                    in
                    if (wsum_new > wmin) then
                        (wsum_new, x::l)
                    else
                        (wsum_new, l)
                  end
              ) (0.,[]) data
        in
        List.rev new_data
    in

    let result =
      match range with
      | (0.,100.)   -> { property ; data }
      | (rmin,rmax) -> { property ; data=data_in_range rmin rmax }
    in
    result






(** Fold function for block values *)
let fold_blocks ~f { property ; data } =
  let init =
    try
      let block = List.hd data in
      Sample.to_float block.Block.value
    with
    | Failure _ -> 0.
  in
  List.fold_left (fun accu block ->
    let x = Sample.to_float block.Block.value
    in f accu x
  ) init data



(** Convergence plot *)
let convergence { property ; data } =

  let rec loop ~sum ~avsq ~ansum ~avsum ~n ~accu = function
    | [] -> List.rev accu
    | head :: tail ->
      begin
        let x = Sample.to_float head.Block.value
        and w = Weight.to_float head.Block.weight
        and avcu0 = avsum /. ansum
        in
        let xw = x *. w
        in
        let ansum = ansum +. w
        and avsum = avsum +. xw
        and sum   = sum   +. xw
        in
        let accu =
          if (n > 0.) then
            (sum /. ansum, sqrt ( abs_float ( avsq /.( ansum *. n))))::accu
          else
            (sum /. ansum, 0.)::accu
        in
        loop tail
          ~sum:sum
          ~avsq:(avsq +. (1. -. (w /. ansum)) *. (x -. avcu0)
                  *. (x -. avcu0) *. w)
          ~avsum:avsum
          ~ansum:ansum
          ~n:(n +. 1.)
          ~accu:accu
      end
  in
  match data with
  | [] -> []
  | head :: tail ->
    begin
      let x = Sample.to_float head.Block.value
      and w = Weight.to_float head.Block.weight
      in
      let s = x *. w in
      loop tail
          ~sum:s
          ~avsq:0.
          ~ansum:w
          ~avsum:s
          ~n:0.
          ~accu:[ (s /. w, 0.) ]
    end


let rev_convergence { property ; data } =
  let p = { property=property ; data = List.rev data } in
  convergence p
  |> List.rev



(** Min and max of block *)
let min_block =
  fold_blocks ~f:(fun accu x ->
    if (x < accu) then x
    else accu
  )


let max_block =
  fold_blocks ~f:(fun accu x ->
    if (x > accu) then x
    else accu
  )



(** Create a hash table for merging *)
let create_hash ~create_key ?(update_block_id=(fun x->x))
    ?(update_value=(fun wc vc wb vb sw -> (wc *. vc +. wb *. vb) /. sw) )
    ?(update_weight=(fun wc wb -> wc +. wb) ) t =
  let table = Hashtbl.create 63
  in
  List.iter (fun block ->
    let key = create_key block
    in
    let open Block in
      hashtbl_change table key (function
      | Some current ->
          let wc, wb =
            Weight.to_float current.weight,
            Weight.to_float block.weight
          in
          let sw =
            update_weight wc wb
          in
          if (Property.is_scalar current.property) then
            let vc, vb =
              Sample.to_float current.value,
              Sample.to_float block.value
            in Some
            {  property = current.property ;
               weight   = Weight.of_float sw ;
               value    = Sample.of_float (update_value wc vc wb vb sw);
               block_id = update_block_id block.block_id;
               pid      = block.pid ;
               compute_node = block.compute_node;
            }
          else
            let vc, vb =
              Sample.to_float_array current.value,
              Sample.to_float_array block.value
            and dim =
              Sample.dimension current.value
            in Some
            {  property = current.property ;
               weight   = Weight.of_float sw ;
               value    =
                 Array.init dim (fun i -> update_value wc vc.(i) wb vb.(i) sw)
                 |> Sample.of_float_array ~dim ;
               block_id = update_block_id block.block_id;
               pid      = block.pid ;
               compute_node = block.compute_node;
            }
      | None -> Some
          { property = block.property ;
            weight   = block.weight;
            value    = block.value ;
            block_id = update_block_id block.block_id;
            pid      = block.pid ;
            compute_node = block.compute_node;
          }
        )
  ) t.data ;
  table



(** Genergic merge function *)
let merge ~create_key ?update_block_id ?update_value ?update_weight t =
  let table = create_hash ~create_key ?update_block_id ?update_value ?update_weight t
  in
  { property = t.property ;
    data = hashtbl_to_alist table
    |>  List.sort (fun x y ->
        if (x>y) then 1
        else if (x<y) then -1
        else 0)
    |>  List.rev_map snd
    |>  List.rev
  }



(** Merge per block id *)
let merge_per_block_id =
  merge
   ~create_key:(fun block -> Block_id.to_string block.Block.block_id)


(** Merge per compute_node *)
let merge_per_compute_node =
  merge
   ~create_key:(fun block ->
      Printf.sprintf "%s"
       (Compute_node.to_string block.Block.compute_node) )



(** Merge per Compute_node and PID *)
let merge_per_compute_node_and_pid =
  merge
   ~create_key:(fun block ->
      Printf.sprintf "%s %10.10d"
       (Compute_node.to_string block.Block.compute_node)
       (block.Block.pid) )



(** Merge per Compute_node and BlockId *)
let merge_per_compute_node_and_block_id =
  merge
   ~create_key:(fun block ->
      Printf.sprintf "%s %10.10d"
       (Compute_node.to_string block.Block.compute_node)
       (Block_id.to_int block.Block.block_id) )

let error_x_over_y = function
  | []          -> (Average.of_float 0., None)
  | (x,_)::[]   -> (Average.of_float x , None)
  | (x,w)::tail ->
    begin
      let avcu0 = ref 0.
      and ansum = ref w
      and avsum = ref x
      and avbl  = ref (x /. w)
      and avsq  = ref 0.
      and n     = ref 1.
      in
      let avcu  = ref !avbl
      in
      List.iter (fun (x,w) ->
          avcu0 := !avsum /. !ansum;
          ansum := !ansum +. w;
          avsum := !avsum +. x;
          avbl  := x /. w ;
          if (!ansum <> 0.) then
            avcu := !avsum /. !ansum
          else ();
          avsq := !avsq +. (1. -. w /. !ansum) *. (!avbl -. !avcu0) *.  (!avbl -. !avcu0) *. w;
          n := !n +. 1.
        ) tail ;
      let arg =
        abs_float (!avsq /.(!ansum *. (!n -. 1.)))
      in
      let error =
        sqrt arg
      in
      (Average.of_float !avcu, Some (Error.of_float error) )
    end


(** Create float, variable operators *)
let one_variable_operator ~update_value p f =
    { p with
      data = List.rev @@ List.rev_map (fun b -> { b with
        Block.value = Sample.of_float (update_value (Sample.to_float b.Block.value) ) }
      ) p.data }

let ( +@ ) p f = one_variable_operator p f
  ~update_value: (fun x -> x +. f )

let ( *@ ) p f = one_variable_operator p f
  ~update_value: (fun x -> x *. f )

let ( -@ ) p f = one_variable_operator p f
  ~update_value: (fun x -> x -. f )

let ( /@ ) p f = one_variable_operator p f
  ~update_value: (fun x -> x /. f )


(** Create two variable operators *)
let two_variable_operator ~update_value p1 p2 =
  merge
  ~update_value
  ~create_key:(fun block ->
      Printf.sprintf "%s %10.10d %10.10d"
       (Compute_node.to_string block.Block.compute_node)
       (Block_id.to_int block.Block.block_id)
       (block.Block.pid) )
  ~update_weight:(fun wc wb -> wc )
  { property = p1.property ;
    data = List.concat [ p1.data ; p2.data ] }

let ( +! ) = two_variable_operator
  ~update_value: (fun wc vc wb vb sw -> (vc +. vb) )

let ( *! ) = two_variable_operator
  ~update_value: (fun wc vc wb vb sw -> (vc *. vb) )

let ( -! ) = two_variable_operator
  ~update_value: (fun wc vc wb vb sw -> (vc -. vb) )

let ( /! ) = two_variable_operator
  ~update_value: (fun wc vc wb vb sw -> (vc /. vb) )




(** Merge two consecutive blocks *)
let compress =
  merge
   ~create_key:(fun block ->
      Printf.sprintf "%s %10.10d %10.10d"
        (Compute_node.to_string block.Block.compute_node) block.Block.pid
        (((Block_id.to_int block.Block.block_id)+1)/2))
   ~update_block_id:(fun block_id ->
      ((Block_id.to_int block_id)+1)/2
      |> Block_id.of_int )




(** Last value on each compute node (for wall_time) *)
let max_value_per_compute_node t =
  let table = Hashtbl.create 63
  in
  let create_key block =
      Printf.sprintf "%s %10.10d"
       (Compute_node.to_string block.Block.compute_node)
       (block.Block.pid)
  in
  List.iter (fun block ->
    let key = create_key block
    in
    let open Block in
    hashtbl_change table key (function
    | Some current ->
        let vc = Sample.to_float current.value
        and vb = Sample.to_float block.value
        in
        if (vc > vb) then
          Some current
        else
          Some block
    | None -> Some block
      )
  ) t.data ;
  { property = t.property ;
    data = hashtbl_to_alist table
    |>  List.sort (fun x y ->
        if (x>y) then 1
        else if (x<y) then -1
        else 0)
    |>  List.rev_map (fun (x,y) -> y)
    |>  List.rev
  }




(** String representation *)
let to_string p =
  match p.property with
  | Property.Cpu   ->  Printf.sprintf "%s" (Time.string_of_sec (sum p))
  | Property.Wall  ->  Printf.sprintf "%s" (Time.string_of_sec (sum (max_value_per_compute_node p)))
  | Property.Accep ->  Printf.sprintf "%16.10f" (average p |> Average.to_float)
  | _ ->
    begin
      if Property.is_scalar p.property then
          match ave_error p with
          | (ave, Some error) ->
              let (ave, error) =
                Average.to_float ave,
                Error.to_float error
              in
              Printf.sprintf "%16.10f +/- %16.10f" ave error
          | (ave, None) ->
              let ave =
                Average.to_float ave
              in
              Printf.sprintf "%16.10f" ave
      else
begin
          match ave_error p with
          | (ave, Some error) ->
              let idxmax =
                Average.dimension ave
              in
              let rec f accu idx =
                if (idx < idxmax) then
                    let (ave, error) =
                      Average.to_float ~idx ave,
                      Error.to_float ~idx error
                    in
                    let s =
                      Printf.sprintf "%8d :  %16.10f +/- %16.10f ;" (idx+1) ave error
                    in
                    (f [@tailcall]) (s :: accu) (idx+1)
                else
                    List.rev (" ]" :: accu)
              in
              f ["[ \n"] 0
              |> String.concat "\n"
          | (ave, None) ->
              Average.to_float ave
              |> Printf.sprintf "%16.10f%!"
end
    end




(** Compress block files : Merge all the blocks computed on the same host *)
let compress_files () =

  Block._raw_data := None;

  let properties =
    Lazy.force Block.properties
  in

  (* Create temporary file *)
  let dir_name =
    Block.dir_name
  in

  let dir_name =
    Lazy.force dir_name
  in
  let files =
    Sys.readdir dir_name
    |> Array.to_list
    |> List.filter (fun x ->
         try
           Str.search_backward (Str.regexp "locked") x (String.length x) >= 0
         with
         | Not_found -> true
       )
    |> List.rev_map (fun x -> dir_name^x)
    |> List.rev
  in

  let out_channel_dir =
    let rand_num = Random.int 1000000 |> string_of_int in
    let dirname =
      Filename.concat !Ezfio.ezfio_filename "blocks"
    in
    if not ( Sys.file_exists dirname ) then
        Unix.mkdir dirname 0o755;
    let tmp_dir =
      Filename.concat dirname ("qmc"^rand_num)
    in
    try
      Unix.mkdir tmp_dir 0o755;
      tmp_dir
    with _ ->
    let message = Printf.sprintf "Cannot create temp dir %s" tmp_dir in
    raise (Sys_error message)
  in

  let out_channel_name =
    let hostname =
      Lazy.force Qmcchem_config.hostname
    and suffix =
     Unix.getpid ()
     |> string_of_int
    in
    String.concat "." [ hostname ; suffix ]
  in

  let block_channel =
    Filename.concat out_channel_dir out_channel_name
    |> open_out
  in

  List.iter (fun p ->
    let l =
      match p with
        | Property.Cpu
        | Property.Accep ->
          of_raw_data ~locked:false ~range:(0.,100.) ~clean:None p
            |> merge_per_compute_node
        | Property.Wall  ->
          of_raw_data ~locked:false ~range:(0.,100.) ~clean:None p
            |> max_value_per_compute_node
        | _     ->
          of_raw_data ~locked:false ~range:(0.,100.) ~clean:None p
            (*
            |> merge_per_compute_node_and_block_id
            *)
    in
    List.iter (fun x ->
      output_string block_channel (Block.to_string_or_bytes x);
      if not Qmcchem_config.binary_io then
        output_char   block_channel '\n';
    ) l.data
  ) properties ;
  close_out block_channel;

  List.iter Unix.unlink files ;
  Unix.rename (Filename.concat out_channel_dir out_channel_name)
              (Filename.concat dir_name        out_channel_name);
  Unix.rmdir  out_channel_dir



(** Autocovariance function (not weighted) *)
let autocovariance { property ; data } =
  let ave =
    average { property ; data }
    |> Average.to_float
  and data =
    match (merge_per_block_id { property ; data })
    with { property ; data } -> Array.of_list data
  in
  let x_t =
     Array.map (fun x -> (Sample.to_float x.Block.value) -. ave) data
  in
  let f i =
    let denom =
      if (i > 1) then (float_of_int i) else 1.
    in
    let r =
      Array.sub x_t 0 i
      |> Array.fold_left (fun accu x ->
        accu +. x *. x_t.(i)) 0.
    in
    r /. denom
  in
  Array.init (Array.length data) f
  |> Array.to_list






(** Computes a histogram *)
let histogram { property ; data } =
  let min, max =
    (min_block { property ; data }),
    (max_block { property ; data })
  in
  let length =
    max -. min
  and n =
    List.length data
    |> float_of_int
    |> sqrt
  in
  let delta_x =
    length /. (n-.1.)
  and result =
    Array.init (int_of_float n + 1) (fun _ -> 0.)
  in
  List.iter (fun x ->
    let w =
      (Weight.to_float x.Block.weight)
    and x =
      (Sample.to_float x.Block.value)
    in
    let i =
      (x -. min) /. delta_x +. 0.5
      |> int_of_float
    in
    result.(i) <- result.(i) +. w
  ) data
  ;
  let norm =
    1. /. ( delta_x *. (
      Array.fold_left (fun accu x -> accu +. x) 0. result
    ) )
  in
  Array.mapi (fun i x -> (min +. (float_of_int i)*.delta_x, x *. norm) ) result
  |> Array.to_list


