let simulation_nucl_fitcusp_factor  = lazy(
  let default =
    1.
  in
  if (Ezfio.has_pseudo_do_pseudo ()) then
    if (Ezfio.get_pseudo_do_pseudo ()) then
      0.
    else
      default
  else
    default
)

let electrons_elec_walk_num         = lazy ( 100         )
let electrons_elec_walk_num_tot     = lazy ( 1000        )
let jastrow_jast_type               = lazy ( "None"      )
let jastrow_jpsi_type               = lazy ( "None"      )
let simulation_block_time           = lazy ( 30          )
let simulation_ci_threshold         = lazy ( 1.e-8       )
let simulation_method               = lazy ( "VMC"       )
let simulation_sampling             = lazy ( "Langevin"  )
let simulation_stop_time            = lazy ( 3600        )
let simulation_time_step            = lazy ( 0.15        )
let simulation_srmc_projection_time = lazy ( 1.          )

let reset_defaults () =
  List.iter (fun x -> Sys.remove ( (Lazy.force Qputils.ezfio_filename) ^ x))
    [ "/electrons/elec_walk_num"       ;
      "/electrons/elec_walk_num_tot"   ;
      "/jastrow/jast_type"             ;
      "/jastrow/jpsi_type"             ;
      "/simulation/block_time"         ;
      "/simulation/ci_threshold"       ;
      "/simulation/method"             ;
      "/simulation/sampling"           ;
      "/simulation/stop_time"          ;
      "/simulation/time_step"          ;
      "/simulation/nucl_fitcusp_factor" ]

