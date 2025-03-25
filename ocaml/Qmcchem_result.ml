open Qptypes

(** Display a table that can be plotted by gnuplot *)
let display_table ~clean ~range ~filter_func property =
  let p =
    Property.of_string property
    |> Random_variable.of_raw_data ~range ~clean ~filter_func
  in
  let  conv = Random_variable.convergence p
  and rconv = Random_variable.rev_convergence p
  and data  = p.Random_variable.data
  in
  let results =
    List.rev_map2 (fun (val1, err1) (val2,err2) -> (val1, err1, val2, err2)) conv rconv
    |> List.rev
  in
  List.iter2 (fun (val1, err1, val2, err2) block ->
     Printf.printf "%10.6f  %10.6f    %10.6f  %10.6f   %10.6f\n"
     val1 err1 val2 err2 (Sample.to_float block.Block.value)
  ) results data


(** Display a convergence plot of the requested property *)
let display_plot ~clean ~range ~filter_func property =
  print_string ("display_plot "^property^".\n")


(** Display a convergence table of the error *)
let display_err_convergence ~clean ~range ~filter_func property =
  let p =
    Property.of_string property
    |> Random_variable.of_raw_data ~range ~clean ~filter_func
  in
  let rec aux n p =
      match Random_variable.ave_error p with
      | (ave, Some error) ->
          let (ave, error) =
            Random_variable.Average.to_float ave,
            Random_variable.Error.to_float error
          in
          Printf.printf "%10d  %16.10f %16.10f\n" n ave error ;
          begin
            if ((3*n) < (List.length p.Random_variable.data)) then
              let new_p =
                 Random_variable.compress p
              in
              aux (n+n) new_p
          end
      | (ave, None) -> ()
  in
  aux 1 p


(** Display the centered cumulants of a property *)
let display_cumulants ~clean ~range ~filter_func property =
  let p =
    Property.of_string property
    |> Random_variable.of_raw_data ~range ~clean ~filter_func
  in
  let cum =
    Random_variable.centered_cumulants p
  in
  Printf.printf "Average     = %16.10f\n" cum.(0);
  Printf.printf "Variance    = %16.10f\n" cum.(1);
  Printf.printf "Centered k3 = %16.10f\n" cum.(2);
  Printf.printf "Centered k4 = %16.10f\n" cum.(3);
  Printf.printf "\n%!";
  let n = 1. /. 12. *. cum.(2) *. cum.(2) +.
          1. /. 48. *. cum.(3) *. cum.(3)
  in
  Printf.printf "Non-gaussianity = %16.10f\n" n



(** Display a table for the autocovariance of the property *)
let display_autocovariance ~clean ~range ~filter_func property =
  let p =
    Property.of_string property
    |> Random_variable.of_raw_data ~range ~clean ~filter_func
  in
  Random_variable.autocovariance p
  |> List.iteri (fun i x ->
      Printf.printf "%10d  %16.10f\n" i x)


(** Display a histogram of the property *)
let display_histogram ~clean ~range ~filter_func property =
  let p =
    Property.of_string property
    |> Random_variable.of_raw_data ~range ~clean ~filter_func
  in
  let histo =
    Random_variable.histogram p

  in
  let g =
    Random_variable.GaussianDist.create
      ~mu:(Random_variable.average p)
      ~sigma2:((Random_variable.centered_cumulants p).(1)
                |> Random_variable.Variance.of_float)
  in
  let g =
    Random_variable.GaussianDist.eval ~g
  in
  List.iter ( fun (x,y) ->
      Printf.printf "%16.10f  %16.10f  %16.10f\n" x y (g ~x)) histo
    (*
  and sigma2 =
    (Random_variable.centered_cumulants p).(1)
  and pi =
    acos(-1.)
  in
  let one_over_2sigma2 =
    1. /. ( 2. *. sigma2 )
  and mu =
    Random_variable.average p
  and norm =
    1. /. (sqrt(sigma2 *. 2.*.pi))
  in
  List.rev_map histo ~f:(fun (x,y) ->
    let g =
      norm *. exp(-.((x-.mu)*.(x-.mu)*.one_over_2sigma2))
    in
    (x,y,g)
  )
  |> List.rev
  |> List.iter ~f:(fun (x,y,g) ->
     Printf.printf "%16.10f  %16.10f  %16.10f\n" x y g)
      *)




(** Display a summary of all the computed quantities *)
let display_summary ~clean ~range ~filter_func =

  let properties =
    Lazy.force Block.properties
  and print_property property =
    Printf.printf "%20s : %!" (Property.to_string property);
    begin
      match property with
      | Accep
      | Wall
      | Cpu -> begin
        let p = Random_variable.of_raw_data ~range ~clean ~filter_func property in
        Printf.printf "%s%!" (Random_variable.to_string p)
        end
      | _ -> begin
        let p = Random_variable.of_raw_data ~range ~clean ~filter_func property in
        Printf.printf "%s%!" (Random_variable.to_string p);
        Printf.printf " (%d)%!" (List.length p.Random_variable.data)
        end
    end;
    Printf.printf "\n%!"
  in
  List.iter print_property properties ;


  let cpu =
    Random_variable.of_raw_data ~range ~clean:None ~filter_func:None Property.Cpu
    |> Random_variable.sum
  and wall =
    Random_variable.of_raw_data ~range ~clean:None ~filter_func:None Property.Wall
    |> Random_variable.max_value_per_compute_node
    |> Random_variable.sum
  and total_weight =
    Random_variable.of_raw_data ~range ~clean ~filter_func Property.E_loc
    |> Random_variable.total_weight
  in

  let speedup =
    cpu /. wall
  in
  Printf.printf "%20s : %10.2f x\n" "Speedup" speedup ;
  Printf.printf "%20s : %20.10e\n" "Total weight" total_weight



let run ?a ?c ?e ?h ?t ?p ?rmin ?rmax ?wmax ?clean ezfio_file =

  Qputils.set_ezfio_filename ezfio_file;

  let rmin =
      match rmin with
      | None -> 0.
      | Some x when (float_of_string x < 0.)   -> failwith "rmin should be >= 0"
      | Some x when (float_of_string x > 100.) -> failwith "rmin should be <= 100"
      | Some x -> float_of_string x
  and rmax =
      match rmax with
      | None -> 100.
      | Some x when (float_of_string x < 0.)   -> failwith "rmax should be >= 0"
      | Some x when (float_of_string x > 100.) -> failwith "rmax should be <= 100"
      | Some x -> float_of_string x
  and wmax =
      match wmax with
      | None -> None
      | Some x when (float_of_string x < 0.)   -> failwith "wmax should be >= 0"
      | Some x when (float_of_string x > 100.) -> failwith "wmax should be <= 100"
      | Some x -> Some (float_of_string x)
  and clean =
      match clean with
      | Some x -> Some (float_of_string x)
      | None -> None
  in


  let range =
    (rmin, rmax)
  in

  let l =
    [ (a, display_autocovariance) ;
      (c, display_cumulants) ;
      (e, display_err_convergence) ;
      (h, display_histogram) ;
      (p, display_plot) ;
      (t, display_table) ;
    ]
  in

  let filter_func: (Block.t list -> Block.t list) option =
    (* Filter out blocks with too large weights *)
    let result =
      match wmax with
      | None -> None
      | Some wmax -> Some (fun data ->
          let max_weight =
            List.fold_left (fun accu x ->
              let w = Weight.to_float x.Block.weight in
              if w > accu then w else accu
              ) 0.  data
          in
          let w = wmax *. 0.01 *. max_weight in
          List.filter (fun x -> Weight.to_float x.Block.weight < w) data
      )
    in result
  in

  List.iter (fun (x,func) ->
      match x with
      | Some property -> func ~clean ~range ~filter_func property
      | None -> ()
    ) l;

  if (List.fold_left (fun accu x ->
      match x with
      | (None, _) -> accu && true
      | (Some _,_) -> false
     ) true l
    ) then
    display_summary ~range ~clean ~filter_func


let command () =
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - QMC=Chem command");
    set_description_doc "Displays the results computed in an EZFIO directory.";

    [ { short='a' ; long="autocovariance" ; opt=Optional ;
        doc="Display the autcovariance function of the property";
        arg=With_arg "<string>" ; };

      { short='c' ; long="centered-cumulants" ; opt=Optional ;
        doc="Print the centered cumulants of a property" ;
        arg=With_arg "<string>"; };

      { short='e' ; long="error" ; opt=Optional ;
        doc="Display the convergence of the error of the property by merging blocks";
        arg=With_arg "<string>"; };

      { short='i' ; long="histogram" ; opt=Optional ;
        doc="Display the histogram of the property blocks" ;
        arg=With_arg "<string>"; };

      { short='p' ; long="plot" ; opt=Optional ;
        doc="Display a convergence plot for a property";
        arg=With_arg "<string>"; };

      { short='m' ; long="rmin" ; opt=Optional ;
        doc="Lower bound of the percentage of the total weight to consider (default 0)" ;
        arg=With_arg "<int>"; };

      { short='n' ; long="rmax" ; opt=Optional ;
        doc="Upper bound of the percentage of the total weight to consider (default 100)" ;
        arg=With_arg "<int>"; };

      { short='x' ; long="clean" ; opt=Optional ;
        doc="Clean values which are beyond x.sigma (default None)" ;
        arg=With_arg "<float>"; };

      { short='t' ; long="table" ; opt=Optional ;
        doc="Print a table for the convergence of a property" ;
        arg=With_arg "<string>"; };

      { short='t' ; long="weight_max" ; opt=Optional ;
        doc="Keep only weights below w/100*max_weight" ;
        arg=With_arg "<string>"; };

      anonymous "EZFIO_DIR" Mandatory "EZFIO directory";
    ]

    |> set_specs ;
  end;

  let a = Command_line.get "autocovariance" in
  let c = Command_line.get "centered-cumulants" in
  let e = Command_line.get "error" in
  let h = Command_line.get "histogram" in
  let t = Command_line.get "table" in
  let p = Command_line.get "plot" in
  let rmin = Command_line.get "rmin" in
  let rmax = Command_line.get "rmax" in
  let wmax = Command_line.get "wmax" in
  let clean = Command_line.get "clean" in

  let ezfio_file =
    match Command_line.anon_args () with
    | ezfio_file :: [] -> ezfio_file
    | _ -> (Command_line.help () ; failwith "Inconsistent command line")
  in
  run ?a ?c ?e ?h ?t ?p ?rmin ?rmax ?wmax ?clean ezfio_file




