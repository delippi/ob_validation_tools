#!/bin/bash

home=`pwd`

#RDASApp="/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/PRs/rrfs-workflow.20260526.dev-sci/sorc/RDASApp"
RDASApp="/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/PRs/rrfs-workflow.20260819.rw/sorc/RDASApp"
rrfswx="/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/PRs/rrfs-workflow.20260819.rw"

echo $RDASApp
cd $RDASApp/ush
source load_rdas.sh
source $rrfswx/ush/fix_rrfs_locations.sh

cd $home

#PDY="20240506"
#cyc="01"
cdate="2024050601"
exp="det_c13_3dvar_jedi_test"
if [[ $# -ge 2 ]]; then
  exp="$1"
  cdate="$2"
fi

sleep 5

cp $FIX_RRFS_LOCATION/gsi/RRFS_CONUS_13km/fv3_grid_spec .

PDY="${cdate:0:8}"
cyc="${cdate:8:2}"
version="v1.2.0"
MATSdir="/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/MATS/ob_validation_tools"
RETROdir="$MATSdir/May2024_retro/${version}/${exp}/rrfs.${PDY}/${cyc}"
grid="./fv3_grid_spec"
lat_tol=0.0005
lon_tol=0.0005
pres_tol=0.1
time_tol=10

# Upper-air conventional
yamls=(
adpupa_airTemperature_120
#adpupa_airTemperature_132
adpupa_specificHumidity_120
#adpupa_specificHumidity_132
adpupa_stationPressure_120
adpupa_winds_220
#adpupa_winds_232

#aircar_airTemperature_133
#aircar_specificHumidity_133
#aircar_winds_233

#aircft_airTemperature_130
#aircft_airTemperature_131
#aircft_airTemperature_134
#aircft_airTemperature_135
##aircft_specificHumidity_134
#aircft_winds_230
#aircft_winds_231
#aircft_winds_234
#aircft_winds_235

#proflr_winds_227
#vadwnd_winds_224
#rassda_airTemperature_126
# Surface conventional
#adpsfc_airTemperature_181  #gsd_terrain_match_surfTobs
#adpsfc_airTemperature_183
#adpsfc_airTemperature_187  #gsd_terrain_match_surfTobs
#adpsfc_specificHumidity_181
#adpsfc_specificHumidity_183
#adpsfc_specificHumidity_187
#adpsfc_stationPressure_181
#adpsfc_stationPressure_187
#adpsfc_winds_281
#adpsfc_winds_284
#adpsfc_winds_287
#msonet_airTemperature_188
#msonet_specificHumidity_188
#msonet_stationPressure_188
#msonet_winds_288

#sfcshp_airTemperature_180
#sfcshp_airTemperature_182
#sfcshp_airTemperature_183
#sfcshp_specificHumidity_180
#sfcshp_specificHumidity_182
#sfcshp_specificHumidity_183
#sfcshp_stationPressure_180
#sfcshp_stationPressure_182
#sfcshp_winds_280
#sfcshp_winds_282
#sfcshp_winds_284

#satwnd_winds_243_56
#satwnd_winds_253_56
#satwnd_winds_243_57
#satwnd_winds_253_57
#satwnd_winds_242_174
#satwnd_winds_250_174
#satwnd_winds_252_174
#satwnd_winds_245_270
#satwnd_winds_246_270
#satwnd_winds_247_270
#satwnd_winds_245_272
#satwnd_winds_246_272
#satwnd_winds_247_272

#nexrad_radialVelocity_999
)

for y in "${yamls[@]}"; do

  platform=$(echo $y | cut -d_ -f1)
  varname=$(echo  $y | cut -d_ -f2)
  otype=$(echo    $y | cut -d_ -f3)

  varname2=$varname
  if [[ $varname == "winds" ]]; then
    varname2="windEastward"
  fi

  # Apply your GSI/JEDI otype offsets
  # original logic preserved
  adj=$otype
  if [[ ${platform} == "satwnd" ]]; then
    stype=$(echo    $y | cut -d_ -f4)
    adj=${otype}_${stype}
    lat_tol=0.05
    lon_tol=0.05
    pres_tol=50
    time_tol=600
  fi
  if [[ ${platform} == "proflr" ]]; then
    pres_tol=2
  fi
  if [[ ${platform} == "vadwnd" ]]; then
    pres_tol=10
  fi
  if [[ ${otype} -eq 281 ]]; then
    pres_tol=5
  fi
  #sfctypes=(187 287 181 281 183 284)
  #if [[ " ${sfctypes[@]} " =~ " ${otype} " ]]; then
  if [[ ${platform} == "adpsfc" ]]; then
    lat_tol=0.001 #0.0005
    lon_tol=0.001 #0.0005
    pres_tol=4 #2
    time_tol=20 #10
  fi
  if [[ ${platform} == "sfcshp" ]]; then
    lat_tol=0.001 #0.0005
    lon_tol=0.001 #0.0005
    pres_tol=4 #2
    time_tol=20 #10
  fi
  if [[ ${platform} == "msonet" ]]; then
    lat_tol=0.001 #0.0005
    lon_tol=0.001 #0.0005
    pres_tol=4 #2
    time_tol=20 #10
  fi
  jedi="${RETROdir}/jdiag_${platform}_${varname}_${adj}.nc"
  gsi="${RETROdir}/gdiag_${platform}_${varname}_${adj}.nc"
  out="hofx_compare_${PDY}${cyc}_${platform}_${varname2}_${adj}_${version}_${exp}"

  echo "Running: $platform  $varname2  $adj"
  echo "  GDIAG = $gsi"
  echo "  JDIAG = $jedi"
  echo "  OUT   = $out"

if [[ ${platform} == "nexrad" ]]; then
python hofx_compare.py -g $gsi -j $jedi --var $varname2 --grid $grid -o $out --lat_tol $lat_tol --lon_tol $lon_tol --height_tol 1.0 --time_tol $time_tol #--debug_oefpc
else
python hofx_compare.py -g $gsi -j $jedi --var $varname2 --grid $grid -o $out --lat_tol $lat_tol --lon_tol $lon_tol --pres_tol $pres_tol --time_tol $time_tol #--debug_oefpc
fi

done
