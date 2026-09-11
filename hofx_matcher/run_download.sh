#!/bin/bash

set -euo pipefail

#cdate=2024050601
#expname="det_c13_3dvar_jedi_baseline005"
user="donald.e.lippi"
workdir=`pwd`

if [[ $# -ge 2 ]]; then
  expname="$1"
  cdate="$2"
fi
pdy="${cdate:0:8}"
cyc="${cdate:8:2}"

#for cyc in $(seq -w 01 23); do
  exp0="May2024_retro/v1.2.0/${expname}/prod/rrfs.${pdy}/${cyc}"
  exp="May2024_retro/v1.2.0/${expname}/rrfs.${pdy}/${cyc}"

  echo "making dir $exp"
  mkdir -p "$workdir/$exp"

  echo "rsync-ing $cyc"
  rsync -av \
    "/lfs/h2/emc/ptmp/${user}/$exp0/"* \
    "$workdir/$exp/." || true

  cd "$workdir/$exp"

  echo "unzipping $cyc"
  shopt -s nullglob
  for f in *.nc4.gz; do
    echo "    Gunzipping $f"
    gunzip -f "$f"
  done
  shopt -u nullglob
#done
