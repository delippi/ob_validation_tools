#!/bin/bash

workdir=`pwd`

gexp="det_c13_3dvar_gsi_simple"
jexp="det_c13_3dvar_jedi_simple"
gexp="det_na3km_3dvar_gsi_simple"
jexp="det_na3km_3dvar_jedi_simple"
cdate=2024050601

cd $workdir
./run_download_ecflow.sh --pair $gexp $jexp $cdate

cd $workdir
./run_convert_diags_ecflow.sh $gexp $cdate

cd $workdir
./run_hofx_compare_ecflow.sh $gexp $jexp $cdate
