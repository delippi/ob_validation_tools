#!/bin/bash

path="/lfs/h2/emc/stmp/donald.e.lippi/May2024_retro/v1.2.0/"
#exp="det_c13_3dvar_jedi_baseline005"
exp="det_c13_3dvar_jedi_sfcCorrected"
exp="det_c13_3dvar_jedi_at2m_identity"

cdate="2024050601"

#gsi="$path/$gsi/$cdate/anal_conv_gsi/OUTPUT*"
gsi="$path/$exp/$cdate/anal_conv_gsi/OUTPUT*"
jedi="$path/$exp/$cdate/anal_conv_jedi/OUTPUT*"

gsi_file=$(ls $gsi | head -n 1)
jedi_file=$(ls $jedi | head -n 1)

#echo $gsi_file
#echo $jedi_file

cp $gsi_file gsi.log
cp $jedi_file jedi.log

python plot_cost.py \
  --jedi "$jedi_file" \
  --gsi "$gsi_file" \
  --cdate "$cdate" \
  --gsi_half \
  -o cost_compare_${exp}.png


