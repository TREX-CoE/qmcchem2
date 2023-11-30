BEGIN_SHELL [ /usr/bin/env python3 ]

data = [ \
("electrons_elec_coord_pool_size"  , "integer"       , "" ),
("electrons_elec_coord_pool"       , "real"          , "(elec_num+1,3,elec_coord_pool_size)" ),
("bi_ortho_mos_mo_l_coef"          , "real"          , "(ao_num,mo_tot_num)" ),
("bi_ortho_mos_coef_psi_right"     , "real"          , ""                    ),
("bi_ortho_mos_use_lr"             , "logical"       , ""                    ),
("bi_ortho_mos_coef_psi_left"      , "real"          , ""                    ),
("electrons_elec_fitcusp_radius"   , "real"          , ""                      ),
("electrons_elec_alpha_num"        , "integer"       , ""                      ),
("electrons_elec_beta_num"         , "integer"       , ""                      ),
("electrons_elec_walk_num"         , "integer"       , ""                      ),
("electrons_elec_walk_num_tot"     , "integer"       , ""                      ),
("ao_basis_ao_num"                 , "integer"       , ""                      ),
("ao_basis_ao_prim_num"            , "integer"       , "(ao_num)"              ),
("ao_basis_ao_nucl"                , "integer"       , "(ao_num)"              ),
("ao_basis_ao_power"               , "integer"       , "(ao_num,3)"            ),
("ao_basis_ao_expo"                , "real"          , "(ao_num,ao_prim_num_max)" ),
("ao_basis_ao_coef"                , "real"          , "(ao_num,ao_prim_num_max)" ),
("ao_two_e_erf_ints_mu_erf"        , "real"          , ""                   ),
("tc_keywords_j1b_type"            , "integer"       , ""                   ),
("tc_keywords_j1b_pen"             , "real"          , "(nucl_num)"         ),
("tc_keywords_j1b_pen_coef"        , "real"          , "(nucl_num)"         ),
("tc_keywords_j1b_coeff"           , "real"          , "(nucl_num)"         ),
("tc_keywords_mu_r_ct"             , "real"          , ""                   ),
("jastrow_inv_sgn_jast"            , "logical"       , ""                   ),
("jastrow_jast_a_up_up"            , "real"          , ""                   ),
("jastrow_jast_a_up_dn"            , "real"          , ""                   ),
("jastrow_jast_b_up_up"            , "real"          , ""                   ),
("jastrow_jast_b_up_dn"            , "real"          , ""                   ),
("jastrow_jast_pen"                , "real"          , "(nucl_num)"          ),
("jastrow_jast_eeN_e_a"            , "real"          , ""                   ),
("jastrow_jast_eeN_e_b"            , "real"          , ""                   ),
("jastrow_jast_eeN_N"              , "real"          , "(nucl_num)"          ),
("jastrow_jast_core_a1"            , "real"          , "(nucl_num)"          ),
("jastrow_jast_core_a2"            , "real"          , "(nucl_num)"          ),
("jastrow_jast_core_b1"            , "real"          , "(nucl_num)"          ),
("jastrow_jast_core_b2"            , "real"          , "(nucl_num)"          ),
("jastrow_jast_qmckl_type_nucl_num", "integer"       , ""                    ),
("jastrow_jast_qmckl_type_nucl_vector", "integer"    , "(nucl_num)"          ),
("jastrow_jast_qmckl_rescale_ee"   , "double precision", ""                  ),
("jastrow_jast_qmckl_rescale_en"   , "double precision", "(jast_qmckl_type_nucl_num)" ),
("jastrow_jast_qmckl_aord_num"   , "integer", "" ),
("jastrow_jast_qmckl_bord_num"   , "integer", "" ),
("jastrow_jast_qmckl_cord_num"   , "integer", "" ),
("jastrow_jast_qmckl_a_vector"   , "double precision", "(jast_qmckl_type_nucl_num*(jast_qmckl_aord_num+1))" ),
("jastrow_jast_qmckl_b_vector"   , "double precision", "(jast_qmckl_bord_num+1)" ),
("jastrow_jast_qmckl_c_vector_size"   , "integer", "" ),
("jastrow_jast_qmckl_c_vector"   , "double precision", "(jast_qmckl_c_vector_size)" ),
("jastrow_jast_type"               , "character*(32)", ""                    ),
("jastrow_jpsi_type"               , "character*(32)", ""                    ),
("simulation_stop_time"            , "integer"       , ""                      ),
("simulation_precision"            , "integer"       , ""                      ),
("simulation_equilibration"        , "logical"       , ""                      ),
("simulation_block_time"           , "integer"       , ""                      ),
("simulation_time_step"            , "real"          , ""                      ),
("simulation_srmc_projection_time" , "real"          , ""                      ),
("simulation_method"               , "character*(32)", ""                      ),
("simulation_nucl_fitcusp_factor"  , "real"          , ""                     ),
("simulation_save_data"            , "logical"       , ""                      ),
("simulation_print_level"          , "integer"       , ""                      ),
("simulation_sampling"             , "character*(32)", ""                      ),
("simulation_ci_threshold"         , "double precision"          , ""                      ),
("simulation_http_server"          , "character*(128)", ""                      ),
("simulation_md5_key"              , "character*(32)" , ""                     ),
("simulation_e_ref"                , "double precision" , ""                   ),
("simulation_e_trial"              , "double precision" , ""                   ),
("simulation_do_run"               , "logical       " , ""                   ),
("simulation_use_trexio"           , "logical"           ,  ""),
("simulation_use_qmckl"            , "logical"           ,  ""),
("trexio_trexio_file"              , "character*(128)" ,  ""),
("pseudo_do_pseudo"                , "logical       " , ""                   ),
("spindeterminants_n_svd_alpha"    , "integer"        , ""                   ),
("spindeterminants_n_svd_beta"     , "integer"        , ""                   ),
("spindeterminants_n_svd_coefs"    , "integer"        , ""                   ),
("spindeterminants_psi_svd_alpha"  , "double precision", "(det_alpha_num,n_svd_alpha,n_states)"),
("spindeterminants_psi_svd_beta"   , "double precision", "(det_beta_num,n_svd_beta,n_states)"),
("spindeterminants_psi_svd_coefs"  , "double precision", "(n_svd_coefs,n_states)"),
("dmc_dress_la"                    , "integer"        , ""                   ),
("dmc_dress_lb"                    , "integer"        , ""                   ),
("dmc_dress_ld"                    , "integer"        , ""                   ),
("dmc_dress_lla"                   , "integer"        , ""                   ),
("dmc_dress_llb"                   , "integer"        , ""                   ),
("dmc_dress_lld"                   , "integer"        , ""                   ),
]

data_no_set = [\
("nuclei_nucl_coord" , "real"    , "(nucl_num,3)" ),
("pseudo_ao_pseudo_grid"   , "double precision" , "(ao_num,pseudo_lmax+pseudo_lmax+1,pseudo_lmax-0+1,nucl_num,pseudo_grid_size)"),
("pseudo_mo_pseudo_grid"   , "double precision" , "(ao_num,pseudo_lmax+pseudo_lmax+1,pseudo_lmax-0+1,nucl_num,pseudo_grid_size)"),
("pseudo_pseudo_dz_k"      , "double precision" , "(nucl_num,pseudo_klocmax)"),
("pseudo_pseudo_dz_kl"     , "double precision" , "(nucl_num,pseudo_kmax,pseudo_lmax+1)"),
("pseudo_pseudo_grid_rmax" , "double precision" , ""),
("pseudo_pseudo_grid_size" , "integer"          , ""),
("pseudo_pseudo_klocmax"   , "integer"          , ""),
("pseudo_pseudo_kmax"      , "integer"          , ""),
("pseudo_pseudo_lmax"      , "integer"          , ""),
("pseudo_pseudo_n_k"       , "integer"          , "(nucl_num,pseudo_klocmax)"),
("pseudo_pseudo_n_kl"      , "integer"          , "(nucl_num,pseudo_kmax,pseudo_lmax+1)"),
("pseudo_pseudo_v_k"       , "double precision" , "(nucl_num,pseudo_klocmax)"),
("pseudo_pseudo_v_kl"      , "double precision" , "(nucl_num,pseudo_kmax,pseudo_lmax+1)"),
("spindeterminants_n_det_alpha"              ,  "integer"           ,  ""),
("spindeterminants_n_det_beta"               ,  "integer"           ,  ""),
("spindeterminants_n_det"                    ,  "integer"           ,  ""),
("spindeterminants_n_int"                    ,  "integer"           ,  ""),
("spindeterminants_bit_kind"                 ,  "integer"           ,  ""),
("spindeterminants_n_states"                 ,  "integer"           ,  ""),
("spindeterminants_psi_det_alpha"            ,  "integer*8"         ,  "(N_int*bit_kind/8,det_alpha_num)"),
("spindeterminants_psi_det_beta"             ,  "integer*8"         ,  "(N_int*bit_kind/8,det_beta_num)"),
("spindeterminants_psi_coef_matrix_rows"     ,  "integer"           ,  "(det_num_input)"),
("spindeterminants_psi_coef_matrix_columns"  ,  "integer"           ,  "(det_num_input)"),
("spindeterminants_psi_coef_matrix_values"   ,  "double precision"  ,  "(det_num_input,N_states)"),
("spindeterminants_psi_left_coef_matrix_values",  "double precision"  ,  "(det_num_input,N_states)"),
]

data_trexio = [\
("mo_basis_mo_num"   , "integer" , ""                   , "mo_num_32"        ),
("mo_basis_mo_coef"  , "real"    , "(ao_num,mo_tot_num)", "mo_coefficient_32"),
("nuclei_nucl_num"   , "integer" , ""                   , "nucleus_num_32"   ),
("nuclei_nucl_charge", "real"    , "(nucl_num)"         , "nucleus_charge_32"),
("mo_basis_mo_coef_aux", "real", "(ao_num,mo_tot_num)", "mo_coefficient_32"),
]

data_trexio_no_fail = [\
("nucl_coord_trexio"   , "real", "(3,nucl_num)"       , "nucleus_coord_32" ),
]


def do_subst(t0,d):
  t = t0
  t = t.replace("$X",d[0])
  t = t.replace("$T",d[1])
  t = t.replace("$D",d[2])
  if len(d) == 4:
    t = t.replace("$Y",d[3])
  if d[1].startswith("character"):
    size = d[1].split("*")[1][1:-1]
    u = "character"
  elif d[1].startswith("double precision"):
    u = d[1].replace(" ","_")
    size = "1"
  elif "*" in d[1]:
    size = "1"
    u = d[1].replace("*","")
  else:
     size = "1"
     u = d[1]
  t = t.replace("$U",u)
  if d[2] == "":
    t = t.replace("$S",size)
  else:
    if size == "1":
      t = t.replace("$S","size(res)")
    else:
      t = t.replace("$S","%s*size(res)"%(size))
  provide = ""
  tmp = d[2].replace('(','').replace(')','')
  for i in "+-*/":
    tmp = tmp.replace(i,',')
  for i in tmp.split(','):
    if ":" in i:
      i = i.split(':')[1]
    try:
     eval(i+"+1")
    except NameError:
     provide += "  PROVIDE "+i+"\n"
  t = t.replace("$P",provide)
  print(t)

t0 = """
subroutine get_$X(res)
  implicit none
  BEGIN_DOC
! Calls EZFIO subroutine to get $X
  END_DOC
  $T                             :: res$D
  integer                        :: ierr, i
  logical                        :: exists
  PROVIDE ezfio_filename
  $P
  if (.not.is_worker) then
    call ezfio_has_$X(exists)
    if (exists) then
      call ezfio_get_$X(res)
      call ezfio_free_$X
    else
      call ezfio_set_$X(res)
    endif
  else
    call zmq_ezfio_has('$X',exists)
    if (exists) then
      call zmq_ezfio_get_$U('$X',res,$S)
    endif
  endif

end
"""


t1 = """
subroutine get_$X(res)
  implicit none
  BEGIN_DOC
! Calls EZFIO subroutine to get $X
  END_DOC
  $T                             :: res$D
  integer                        :: ierr
  PROVIDE ezfio_filename
  $P
  if (.not.is_worker) then
    call ezfio_get_$X(res)
    call ezfio_free_$X
  else
    call zmq_ezfio_get_$U('$X',res,$S)
  endif

end
"""

for i,d in enumerate(data):
  do_subst(t0,d)

for i,d in enumerate(data_no_set):
  i += len(data)
  do_subst(t1,d)

t2 = """
subroutine get_$X(res)
  use trexio
  implicit none
  BEGIN_DOC
! Calls EZFIO or TREXIO subroutine to get $X
  END_DOC
  $T                             :: res$D
  integer                        :: ierr
  PROVIDE ezfio_filename
  $P
  if (.not.is_worker) then
    if (use_trexio) then
      ierr = trexio_read_$Y(trexio_file, res)
      if (ierr /= TREXIO_SUCCESS) then
        print *, 'Error in TREXIO. Using EZFIO instead'
        call ezfio_get_$X(res)
        call ezfio_free_$X
      end if
    else
      call ezfio_get_$X(res)
      call ezfio_free_$X
    endif
  else
    call zmq_ezfio_get_$U('$X',res,$S)
  endif

end
"""

for i,d in enumerate(data_trexio):
  do_subst(t2,d)

t3 = """
subroutine get_$X(res)
  use trexio
  implicit none
  BEGIN_DOC
! Calls EZFIO or TREXIO subroutine to get $X
  END_DOC
  $T                             :: res$D
  integer                        :: ierr
  PROVIDE ezfio_filename
  $P
  if (.not.is_worker) then
    if (use_trexio) then
      ierr = trexio_read_$Y(trexio_file, res)
      if (ierr /= TREXIO_SUCCESS) then
        print *, irp_here
        character*(128) :: msg
        call trexio_string_of_error(ierr, msg)
        print *, trim(msg)
        stop -1
      end if
    else
      print *, 'Error: use_trexio should be true'
      stop -1
    endif
  else
    call zmq_ezfio_get_$U('$X',res,$S)
  endif

end
"""

for i,d in enumerate(data_trexio_no_fail):
  do_subst(t3,d)

END_SHELL

