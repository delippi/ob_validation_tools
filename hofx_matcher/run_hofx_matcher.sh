#!/bin/bash

workdir=`pwd`

expname="det_c13_3dvar_jedi_sfcCorrected"
expname="det_c13_3dvar_jedi_at2m_identity"
#expname="det_c13_3dvar_jedi_test"
cdate=2024050601

cd $workdir
./run_download.sh $expname $cdate

cd $workdir
./run_convert_diags.sh $expname $cdate

cd $workdir
if [[ $expname == "det_c13_3dvar_jedi_sfcCorrected" ]]; then
  ./run_sfc2.sh $expname $cdate
elif [[ $expname == "det_c13_3dvar_jedi_at2m_identity" ]]; then
  ./run_sfc3.sh $expname $cdate
else
  ./run_hofx_compare.sh $expname $cdate
fi

mkdir -p hofx_match_$expname
mv *"$expname"*.png hofx_match_$expname/.


cp -r May2024_retro/v1.2.0/$expname hofx_match_$expname/.
