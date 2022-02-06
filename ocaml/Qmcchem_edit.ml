let file_header filename = Printf.sprintf
"
+----------------------------------------------------------------+
|                            QMC=Chem                            |
+----------------------------------------------------------------+

Editing file `%s`

" filename

let make_header s =
  let l = String.length s in
  "\n\n"^s^"\n"^(String.init l (fun _ -> '='))^"\n\n"


type field = 
 | Block_time
 | Walk_num
 | Walk_num_tot
 | Stop_time
 | Fitcusp_factor
 | Method
 | Sampling
 | Ref_energy
 | Trial_wf_energy
 | CI_threshold
 | Time_step
 | SRMC_projection_time 
 | Jastrow_type
 | Properties


let get field = 
  let option_to_string read to_string doc = 
    let value =
      read () |> to_string
    in
    Printf.sprintf "%s ::\n\n    %s\n\n" doc value
  in
  let option_to_string_prop read to_string doc = 
    let value =
      read () |> to_string
    in
    Printf.sprintf "%s :\n\n%s\n\n" doc value
  in
  let open Input in
  match field with
 | Block_time   ->
   option_to_string Block_time.read   Block_time.to_string   Block_time.doc  
 | Walk_num     ->
   option_to_string Walk_num.read     Walk_num.to_string     Walk_num.doc  
 | Walk_num_tot ->
   option_to_string Walk_num_tot.read Walk_num_tot.to_string Walk_num_tot.doc  
 | Stop_time    ->
   option_to_string Stop_time.read    Stop_time.to_string    Stop_time.doc  
 | Fitcusp_factor ->
   option_to_string Fitcusp_factor.read      Fitcusp_factor.to_string      Fitcusp_factor.doc
 | Method       ->
   option_to_string Method.read       Method.to_string       Method.doc  
 | Sampling     ->
   option_to_string Sampling.read     Sampling.to_string     Sampling.doc  
 | Ref_energy   ->
   option_to_string Ref_energy.read   Ref_energy.to_string   Ref_energy.doc  
 | Trial_wf_energy ->
   option_to_string Trial_wf_energy.read   Trial_wf_energy.to_string   Trial_wf_energy.doc  
 | CI_threshold ->                    
   option_to_string CI_threshold.read   CI_threshold.to_string   CI_threshold.doc  
 | Time_step    ->                    
   option_to_string Time_step.read    Time_step.to_string    Time_step.doc  
 | SRMC_projection_time ->             
   option_to_string SRMC_projection_time.read SRMC_projection_time.to_string SRMC_projection_time.doc 
 | Jastrow_type ->                    
   option_to_string Jastrow_type.read Jastrow_type.to_string Jastrow_type.doc  
 | Properties   ->                    
   option_to_string_prop Properties.read Properties.to_string Properties.doc  

                                      

let create_temp_file ?temp_filename ezfio_filename fields =
  let filename =
    match temp_filename with
    | None -> Filename.temp_file "qmcchem_edit_" ".rst"
    | Some name -> name
  in
  let out_channel = open_out filename in
  (file_header ezfio_filename) :: (List.rev @@ List.rev_map get fields)
  |> String.concat "\n"
  |> output_string out_channel
  ; close_out out_channel
  ; filename


(** Write the input file corresponding to the MD5 key *)
let write_input_in_ezfio  ezfio_filename  fields =
  let dirname = 
    Lazy.force  QmcMd5.input_directory
  in
  let temp_filename = 
    QmcMd5.hash ()
    |> Filename.concat  dirname
  in
  let input_filename =
    create_temp_file  ~temp_filename  ezfio_filename  fields
  in
  assert (Sys.file_exists  input_filename)


(** Run the edit command *)
let run ~c ?f ?t ?l ?m ?e ?et ?s ?ts ?w ?wt ?n ?j ?p ?input ezfio_filename =

  let interactive = ref (
    if c then
      false
    else
      true
  )
  in

  (* Open EZFIO *)
  Qputils.set_ezfio_filename ezfio_filename;

  let handle_option (type_conv, write) x =
    let () =
      match x with
      | Some x -> 
        begin
          type_conv x |> write;
          interactive := false;
        end
      | None   -> ()
    in ();
  in

  handle_option      Input.Ref_energy.(of_string, write) e;
  handle_option Input.Trial_wf_energy.(of_string, write) et;
  handle_option    Input.Jastrow_type.(of_string, write) j;
  handle_option      Input.Block_time.(of_string, write) l;
  handle_option          Input.Method.(of_string, write) m;
  handle_option       Input.Stop_time.(of_string, write) t;
  handle_option        Input.Sampling.(of_string, write) s;
  handle_option  Input.Fitcusp_factor.(of_string, write) f;
  handle_option       Input.Time_step.(of_string, write) ts;
  handle_option        Input.Walk_num.(of_string, write) w;
  handle_option    Input.Walk_num_tot.(of_string, write) wt;
  handle_option    Input.CI_threshold.(of_string, write) n;
  handle_option Input.SRMC_projection_time.(of_string, write) p;
                 

  let fields = 
   [ 
      Stop_time            ;
      Block_time           ;
      Method               ;
      Sampling             ;
      Time_step            ;
      SRMC_projection_time  ;
      Ref_energy           ;
      Trial_wf_energy      ;
      Walk_num             ;
      Walk_num_tot         ;
      Fitcusp_factor       ;
      CI_threshold         ;
      Jastrow_type         ;
      Properties           ;
   ]
  in

  if (!interactive) then
    begin
      let temp_filename = 
          create_temp_file ezfio_filename fields
      in
      let () =
        match input with
        | Some filename ->
          begin
            if (not !interactive) then
               failwith "Input file not allowed with command line arguments"
            else
              let rc = 
                Printf.sprintf "cp %s %s" filename temp_filename
                |> Sys.command 
              in
              assert (rc = 0)
          end
        | None          -> 
          begin
            (* Open the temp file with external editor *)
            let editor =
              try Sys.getenv "EDITOR" with
              | Not_found -> "vi"
            in
            let rc = 
              Printf.sprintf "%s %s ; tput sgr0 2> /dev/null" editor temp_filename 
              |> Sys.command
            in
            assert (rc = 0)
          end
      in
     
      (* Re-read the temp file *)
      let re_data =
         Str.regexp "   .+ *$"
      and re_prop =
         Str.regexp "([ xX]) .*$"
      and raw_data =
        let ic = open_in temp_filename in
        let result = String_ext.input_lines ic in
        close_in ic ; result
      in
      let data = 
        ( List.filter (fun x -> Str.string_match re_data x 0) raw_data
          |> List.rev_map String.trim |> List.rev ) @
        [
        List.filter (fun x -> Str.string_match re_prop x 0) raw_data 
        |> List.rev_map String.trim
        |> List.rev
        |> String.concat "\n" ]
      in
      let open Input in
      List.iter2 (fun s f -> 
       try 
         begin
            match f with
            |  Stop_time             ->  Stop_time.(of_string            s  |>  write)
            |  Fitcusp_factor        ->  Fitcusp_factor.(of_string       s  |>  write)
            |  Block_time            ->  Block_time.(of_string           s  |>  write)
            |  Method                ->  Method.(of_string               s  |>  write)
            |  Ref_energy            ->  Ref_energy.(of_string           s  |>  write)
            |  Sampling              ->  Sampling.(of_string             s  |>  write)
            |  Time_step             ->  Time_step.(of_string            s  |>  write)
            |  SRMC_projection_time  ->  SRMC_projection_time.(of_string s  |>  write)
            |  Walk_num              ->  Walk_num.(of_string             s  |>  write)
            |  Walk_num_tot          ->  Walk_num_tot.(of_string         s  |>  write)
            |  CI_threshold          ->  CI_threshold.(of_string         s  |>  write)
            |  Jastrow_type          ->  Jastrow_type.(of_string         s  |>  write)
            |  Trial_wf_energy       ->  Trial_wf_energy.(of_string      s  |>  write)
            |  Properties            ->  Properties.(of_string           s  |>  write)
         end
       with
       | Failure msg -> Printf.eprintf "%s\n" msg
      ) data fields ;
     
      (* Remove temp_file *)
      Sys.remove temp_filename;

    end
  ;

  if c then
    begin
      let dirname = 
        Filename.concat (Filename.concat ezfio_filename "blocks") (QmcMd5.hash ()) 
      in
      let rec clean_dir y =
        if Sys.file_exists y && Sys.is_directory y then
          begin
            Sys.readdir y
            |> Array.map (fun x -> Filename.concat y x)
            |> Array.iter (function x ->
                if Sys.file_exists x && Sys.is_directory x  then
                  clean_dir x
                else
                  Sys.remove x
              );
            Unix.rmdir y
          end
      in clean_dir dirname;
      Printf.printf "Blocks cleared\n"
    end
  ;

  Input.validate ();
  QmcMd5.reset_hash ();
  write_input_in_ezfio ezfio_filename fields


let command () = 
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - QMC=Chem command");
    set_description_doc "Edits input data";

    [ { short='c' ; long="clear" ; opt=Optional ;
        doc="Clears blocks" ;
        arg=Without_arg ; };

      { short='e' ; long="ref-energy" ; opt=Optional ;
        doc=Input.Ref_energy.doc;
        arg=With_arg "<float>"; };
    
      { short='f' ; long="fitcusp" ; opt=Optional ;
        doc=Input.Fitcusp_factor.doc;
        arg=With_arg "<float>"; };
    
      { short='i' ; long="time-step" ; opt=Optional ;
        doc=Input.Time_step.doc;
        arg=With_arg "<float>"; };
    
      { short='j' ; long="jastrow" ; opt=Optional ;
        doc=Input.Jastrow_type.doc;
        arg=With_arg "<string>"; };
    
      { short='l' ; long="block-time" ; opt=Optional ;
        doc=Input.Block_time.doc;
        arg=With_arg "<int>"; };
    
      { short='m' ; long="method" ; opt=Optional ;
        doc=Input.Method.doc;
        arg=With_arg "<string>"; };
    
      { short='n' ; long="norm" ; opt=Optional ;
        doc=Input.CI_threshold.doc;
        arg=With_arg "<float>"; };
    
      { short='p' ; long="projection-time" ; opt=Optional ;
        doc=Input.SRMC_projection_time.doc;
        arg=With_arg "<float>"; };
    
      { short='r' ; long="trial-energy" ; opt=Optional ;
        doc=Input.Trial_wf_energy.doc;
        arg=With_arg "<float>"; };
    
      { short='s' ; long="sampling" ; opt=Optional ;
        doc=Input.Sampling.doc;
        arg=With_arg "<string>"; };
    
      { short='t' ; long="stop-time" ; opt=Optional ;
        doc=Input.Stop_time.doc;
        arg=With_arg "<int>"; };
    
      { short='w' ; long="walk-num" ; opt=Optional ;
        doc=Input.Walk_num.doc;
        arg=With_arg "<int>"; };
    
      { short='x' ; long="walk-num-tot" ; opt=Optional ;
        doc=Input.Walk_num_tot.doc;
        arg=With_arg "<int>"; };
    
      anonymous "EZFIO_DIR" Mandatory "EZFIO directory";
      anonymous "FILE"   Optional  "Name of the input file";
    ]
    |> set_specs 
  end;        

  let c = Command_line.get_bool "clear" in
  let f = Command_line.get "fitcusp" in
  let t = Command_line.get "stop-time" in
  let l = Command_line.get "block-time" in
  let m = Command_line.get "method" in
  let e = Command_line.get "ref-energy" in
  let et = Command_line.get "trial-energy" in
  let s = Command_line.get "sampling" in
  let ts = Command_line.get "time-step" in
  let w = Command_line.get "walk-num" in
  let wt = Command_line.get "walk-num-tot" in
  let n = Command_line.get "norm" in
  let j = Command_line.get "jastrow" in
  let p = Command_line.get "projection-time" in

  let ezfio_file, input =
    match Command_line.anon_args () with
    | ezfio_file :: [] -> ezfio_file, None
    | ezfio_file :: file :: [] -> ezfio_file, (Some file)
    | _ -> (Command_line.help () ; failwith "Inconsistent command line")
  in
  run ~c ?f ?t ?l ?m ?e ?et ?s ?ts ?w ?wt ?n ?j ?p ?input ezfio_file




