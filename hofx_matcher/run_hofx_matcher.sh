#!/bin/bash

workdir=`pwd`

expname="det_c13_3dvar_jedi_baseline005"
cdate=2024050601

cd $workdir
#./run_download.sh $expname $cdate

cd $workdir
#./run_convert_diags.sh $expname $cdate

cd $workdir
./run_hofx_compare.sh $expname $cdate
