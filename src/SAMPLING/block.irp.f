! Providers of *_block_walk
!==============================
BEGIN_SHELL [ /usr/bin/env python3 ]
from properties import *

t = """
 BEGIN_PROVIDER [ $T, $X_block_walk $D1 ]
&BEGIN_PROVIDER [ $T, $X_2_block_walk $D1 ]
 implicit none
 include '../types.F'
 BEGIN_DOC
! Properties averaged over the block per walker
 END_DOC

 if (qmc_method == t_VMC) then
  PROVIDE E_loc_vmc_block_walk
  if (calc_$X) then
    $X_block_walk = $X_vmc_block_walk
    $X_2_block_walk = $X_2_vmc_block_walk
  endif
 else if (qmc_method == t_DMC) then
  PROVIDE E_loc_dmc_block_walk
  if (calc_$X) then
    $X_block_walk = $X_dmc_block_walk
    $X_2_block_walk = $X_2_dmc_block_walk
  endif
 else if (qmc_method == t_SRMC) then
  PROVIDE E_loc_srmc_block_walk
  if (calc_$X) then
    $X_block_walk = $X_srmc_block_walk
    $X_2_block_walk = $X_2_srmc_block_walk
  endif
 else if (qmc_method == t_PDMC) then
  PROVIDE E_loc_pdmc_block_walk
  if (calc_$X) then
    $X_block_walk = $X_pdmc_block_walk
    $X_2_block_walk = $X_2_pdmc_block_walk
  endif
 else if (qmc_method == t_FKMC) then
  PROVIDE E_loc_fkmc_block_walk
  if (calc_$X) then
    $X_block_walk = $X_fkmc_block_walk
    $X_2_block_walk = $X_2_fkmc_block_walk
  endif
 endif

END_PROVIDER
"""

for p in properties:
    if p[2] == "":
      D1 = ""
    else:
      D1 = ", ("+p[2][1:-1]+")"
    print(t.replace("$X",p[1]).replace("$T",p[0]).replace("$D1",D1))
END_SHELL




BEGIN_PROVIDER [ double precision, block_weight ]
 implicit none
 include '../types.F'
 BEGIN_DOC  
! Weight of the current block in the full average of the simulation
 END_DOC
 integer :: i
 block_weight = 0.d0
END_PROVIDER


