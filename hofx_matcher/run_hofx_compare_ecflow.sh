#!/bin/bash

set -eo pipefail
shopt -s nullglob

# Uncomment the diagnostic keys to compare. Leave the list empty to compare
# every key present in both jobs. Keys omit the gdiag_/jdiag_ prefix and .nc.
yamls=(
  adpupa_airTemperature_120
#  #adpupa_airTemperature_132
  adpupa_specificHumidity_120
#  #adpupa_specificHumidity_132
  adpupa_stationPressure_120
  adpupa_winds_220
#  #adpupa_winds_232
#
  aircar_airTemperature_133
  aircar_specificHumidity_133
  aircar_winds_233
#
#  #aircft_airTemperature_130
#  #aircft_airTemperature_131
#  #aircft_airTemperature_134
#  #aircft_airTemperature_135
#  #aircft_specificHumidity_134
#  #aircft_winds_230
#  #aircft_winds_231
#  #aircft_winds_234
#  #aircft_winds_235
#
#  #proflr_winds_227
#  #vadwnd_winds_224
#  #rassda_airTemperature_126

#  # Surface conventional
  adpsfc_airTemperatureAt2M_181  #gsd_terrain_match_surfTobs
  adpsfc_airTemperatureAt2M_183
  adpsfc_airTemperatureAt2M_187  #gsd_terrain_match_surfTobs
  msonet_airTemperatureAt2M_188
  adpsfc_specificHumidityAt2M_181
  adpsfc_specificHumidityAt2M_183
  adpsfc_specificHumidityAt2M_187
  msonet_specificHumidityAt2M_188
  adpsfc_stationPressureAt2M_181
  adpsfc_stationPressureAt2M_187
  msonet_stationPressureAt2M_188
  adpsfc_winds_281
  adpsfc_winds_287
  msonet_winds_288
#
#  #sfcshp_specificHumidityAt2M_180
#  #sfcshp_specificHumidityAt2M_182
#  #sfcshp_specificHumidityAt2M_183
#  #sfcshp_stationPressureAt2M_180
#  #sfcshp_stationPressureAt2M_182
#  #sfcshp_winds_280
#  #sfcshp_winds_282
#  #sfcshp_winds_284
#
#  #satwnd_winds_243_56
#  #satwnd_winds_253_56
#  #satwnd_winds_243_57
#  #satwnd_winds_253_57
#  #satwnd_winds_242_174
#  #satwnd_winds_250_174
#  #satwnd_winds_252_174
#  #satwnd_winds_245_270
#  #satwnd_winds_246_270
#  #satwnd_winds_247_270
#  #satwnd_winds_245_272
#  #satwnd_winds_246_272
#  #satwnd_winds_247_272
#
#  #nexrad_radialVelocity_999
)

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 GSI_EXP_NAME JEDI_EXP_NAME YYYYMMDDHH" >&2
  exit 2
fi

gexp=$1
jexp=$2
cdate=$3
exp_pattern='^[[:alnum:]_][[:alnum:]_.-]*$'
if [[ ! $gexp =~ $exp_pattern || ! $jexp =~ $exp_pattern ||
      ! $cdate =~ ^[0-9]{8}([01][0-9]|2[0-3])$ ]]; then
  echo "Invalid experiment name or cycle." >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
workdir=$PWD
pdy=${cdate:0:8}
cyc=${cdate:8:2}

# Print the selected path on stdout and keep the interactive menu on stderr.
select_job() {
  local exp=$1 task=$2 candidate run_id
  local cycle_dir="ecflow/$exp/rrfs.$pdy/$cyc"
  local run_prefix="rrfs_det_${task}_$cyc"
  local -a run_ids=()
  local PS3='Select a job number: '

  for candidate in "$cycle_dir/$run_prefix".*; do
    run_id=${candidate##*/$run_prefix.}
    if [[ -d $candidate && $run_id =~ ^[0-9]+\.[[:alnum:]_][[:alnum:]_.-]*$ ]]; then
      run_ids+=("$run_id")
    fi
  done
  if (( ${#run_ids[@]} == 0 )); then
    echo "No staged $task jobs in: $cycle_dir" >&2
    return 1
  fi

  mapfile -t run_ids < <(printf '%s\n' "${run_ids[@]}" | LC_ALL=C sort -t . -k1,1n)
  printf 'Staged %s jobs for %s %s\n' "$task" "$exp" "$cdate" >&2
  select run_id in "${run_ids[@]}"; do
    if [[ -n $run_id ]]; then
      printf '%s/%s.%s\n' "$cycle_dir" "$run_prefix" "$run_id"
      return 0
    fi
    echo "Choose a job number from the list." >&2
  done
  echo "No $task job selected." >&2
  return 1
}

gsi_dir=$(select_job "$gexp" analysis_gsi_diag)
jedi_dir=$(select_job "$jexp" analysis_jedi)

# Pair converted GSI and native JEDI diagnostics by their common filename key.
gsi_diags=("$gsi_dir"/gdiag_*.nc)
if (( ${#gsi_diags[@]} == 0 )); then
  echo "No converted GSI diagnostics in $gsi_dir. Run run_convert_diags_ecflow.sh for this job first." >&2
  exit 1
fi

keys=("${yamls[@]}")
if (( ${#keys[@]} == 0 )); then
  for gdiag in "${gsi_diags[@]}"; do
    key=${gdiag##*/gdiag_}
    keys+=("${key%.nc}")
  done
fi

pairs=()
for key in "${keys[@]}"; do
  if [[ -f $gsi_dir/gdiag_$key.nc && -f $jedi_dir/jdiag_$key.nc ]]; then
    pairs+=("$key")
  elif (( ${#yamls[@]} > 0 )); then
    printf 'Skipping %s: missing a GSI or JEDI diagnostic.\n' "$key" >&2
  fi
done
if (( ${#pairs[@]} == 0 )); then
  echo "No selected diagnostic pairs found between $gsi_dir and $jedi_dir." >&2
  exit 1
fi

# Reuse the existing comparison environment on WCOSS2 when needed.
if [[ -n ${RDASApp:-} ]] ||
   ! python -c 'import cartopy, matplotlib, netCDF4, numpy, scipy' >/dev/null 2>&1; then
  RDASApp=${RDASApp:-/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/PRs/rrfs-workflow.20260819.rw/sorc/RDASApp}
  if [[ ! -r "$RDASApp/ush/load_rdas.sh" ]]; then
    echo "Set RDASApp to your RDASApp build or activate Python with cartopy, matplotlib, netCDF4, numpy, and scipy." >&2
    exit 1
  fi
  cd -- "$RDASApp/ush"
  source ./load_rdas.sh
  cd -- "$workdir"
fi

# The grid is optional and only supplies the plotted domain outline.
grid=${HOFX_GRID:-}
if [[ -z $grid && -f fv3_grid_spec ]]; then
  grid=fv3_grid_spec
elif [[ -z $grid && -n ${FIX_RRFS_LOCATION:-} &&
        -f $FIX_RRFS_LOCATION/gsi/RRFS_CONUS_13km/fv3_grid_spec ]]; then
  grid="$FIX_RRFS_LOCATION/gsi/RRFS_CONUS_13km/fv3_grid_spec"
fi
grid_args=()
if [[ -n $grid ]]; then
  if [[ ! -r $grid ]]; then
    echo "Cannot read grid file: $grid" >&2
    exit 1
  fi
  grid_args=(--grid "$grid")
fi

# Keep plots from different experiment and job combinations separate.
out_dir="hofx_match_ecflow/$cdate/$gexp/${gsi_dir##*/}/$jexp/${jedi_dir##*/}"
mkdir -p -- "$out_dir"
plot_tmp=$(mktemp -d "$workdir/$out_dir/.plots.XXXXXX")
trap 'rm -rf -- "$plot_tmp"' EXIT
plots=0
failed=0
for key in "${pairs[@]}"; do
  IFS=_ read -r platform varname otype <<< "$key"
  if [[ $varname == winds ]]; then
    varname=windEastward
  fi

  # Apply the legacy matching tolerances independently to each diagnostic.
  lat_tol=0.0005
  lon_tol=0.0005
  pres_tol=0.1
  time_tol=10
  case "$platform" in
    satwnd)
      lat_tol=0.05
      lon_tol=0.05
      pres_tol=50
      time_tol=600
      ;;
    proflr) pres_tol=2 ;;
    vadwnd) pres_tol=10 ;;
    adpsfc|sfcshp|msonet)
      lat_tol=0.001
      lon_tol=0.001
      pres_tol=4
      time_tol=20
      ;;
  esac
  if [[ ${otype%%_*} == 281 && $platform != adpsfc &&
        $platform != sfcshp && $platform != msonet ]]; then
    pres_tol=5
  fi

  out="hofx_compare_${cdate}_$key"
  printf 'Comparing %s\n' "$key"
  if ! python "$script_dir/hofx_compare.py" \
    -g "$gsi_dir/gdiag_$key.nc" -j "$jedi_dir/jdiag_$key.nc" \
    --var "$varname" "${grid_args[@]}" -o "$plot_tmp/$out" \
    --lat_tol "$lat_tol" --lon_tol "$lon_tol" --pres_tol "$pres_tol" \
    --height_tol 1.0 --time_tol "$time_tol"; then
    echo "Comparison failed for $key." >&2
    failed=$((failed + 1))
    continue
  fi

  # The Python tool exits successfully without a plot when no obs match.
  if [[ -f $plot_tmp/${out}_8panel.png ]]; then
    mv -f -- "$plot_tmp/${out}_8panel.png" "$out_dir/"
    plots=$((plots + 1))
  fi
done

printf 'Wrote %s plot(s) to %s\n' "$plots" "$out_dir"
if (( failed > 0 || plots == 0 )); then
  exit 1
fi

