#!/bin/bash

path="/lfs/h2/emc/stmp/donald.e.lippi/May2024_retro/v1.2.0/"
#gsi="det_c13_3dvar_gsi_t187"
exp="det_c13_3dvar_jedi_test"
#gsi="det_c13_3dvar_gsi_car_qopt2_2_nobal"
#jedi="det_c13_3dvar_gsi_car_qopt2_2_nobal"

cdate="2024050606"
cdate="2024051018"
cdate="2024050601"

#gsi="$path/$gsi/$cdate/anal_conv_gsi/OUTPUT*"
log1="$path/$exp/$cdate/ninner50/OUTPUT*"
log2="$path/$exp/$cdate/ninner30/OUTPUT*"

log1_file=$(ls $log1 | head -n 1)
log2_file=$(ls $log2 | head -n 1)

python plot_cost_same.py \
  --log1 $log1_file \
  --log2 $log2_file \
  --label1 "ninner=50" \
  --label2 "ninner=30" \
  --cdate "$cdate" \
  -o cost_compare_same.png
