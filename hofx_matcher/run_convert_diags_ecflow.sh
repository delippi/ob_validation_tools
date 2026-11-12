#!/bin/bash

set -eo pipefail
shopt -s nullglob

run_id_pattern='[0-9]+\.[[:alnum:]_][[:alnum:]_.-]*'

case $# in
  1)
    run_dir=$1
    ;;
  2)
    expname=$1
    cdate=$2
    if [[ ! $expname =~ ^[[:alnum:]_][[:alnum:]_.-]*$ ||
          ! $cdate =~ ^[0-9]{8}([01][0-9]|2[0-3])$ ]]; then
      echo "Invalid experiment name or cycle." >&2
      exit 2
    fi

    # Match the downloader's staging layout and offer only GSI diagnostic jobs.
    cyc=${cdate:8:2}
    cycle_dir="ecflow/$expname/rrfs.${cdate:0:8}/$cyc"
    run_prefix="rrfs_det_analysis_gsi_diag_$cyc"
    run_ids=()
    for candidate in "$cycle_dir/$run_prefix".*; do
      run_id=${candidate##*/$run_prefix.}
      if [[ -d $candidate && $run_id =~ ^${run_id_pattern}$ ]]; then
        run_ids+=("$run_id")
      fi
    done

    if (( ${#run_ids[@]} == 0 )); then
      echo "No staged GSI diagnostic jobs in: $cycle_dir" >&2
      exit 1
    fi

    # Sort by numeric job ID, oldest first, just as in the downloader.
    mapfile -t run_ids < <(printf '%s\n' "${run_ids[@]}" | LC_ALL=C sort -t . -k1,1n)
    printf 'Staged GSI diagnostic jobs for %s %s\n' "$expname" "$cdate" >&2
    PS3='Select a GSI diagnostic job number: '
    run_id=
    select run_id in "${run_ids[@]}"; do
      if [[ -n $run_id ]]; then
        break
      fi
      echo "Choose a job number from the list." >&2
    done

    if [[ -z $run_id ]]; then
      echo "No job selected." >&2
      exit 1
    fi
    run_dir="$cycle_dir/$run_prefix.$run_id"
    ;;
  *)
    printf 'Usage: %s GSI_EXP_NAME YYYYMMDDHH\n' "$0" >&2
    printf '       %s STAGED_GSI_DIAG_JOB_DIR\n' "$0" >&2
    exit 2
    ;;
esac

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$run_dir"

# Recover the cycle from the directory layout created by the ecFlow downloader.
run_name=${PWD##*/}
cycle_dir=${PWD%/*}
cyc=${cycle_dir##*/}
day_dir=${cycle_dir%/*}
day_dir=${day_dir##*/}
if [[ ! $day_dir =~ ^rrfs\.[0-9]{8}$ ||
      ! $cyc =~ ^([01][0-9]|2[0-3])$ ||
      ! $run_name =~ ^rrfs_det_analysis_gsi_diag_${cyc}\.${run_id_pattern}$ ]]; then
  echo "Expected a staged GSI diagnostic job directory: $PWD" >&2
  exit 2
fi
cdate="${day_dir#rrfs.}$cyc"

# The Python converter takes analysis files and opens the corresponding guess
# files itself. Only require a pair for variables present in this job.
analysis_files=()
for var in t q ps uv rw; do
  anl_file="diag_conv_${var}_anl.${cdate}.nc4"
  ges_file="diag_conv_${var}_ges.${cdate}.nc4"
  if [[ -f $anl_file ]]; then
    if [[ ! -f $ges_file ]]; then
      echo "Missing guess diagnostic for $anl_file: $PWD/$ges_file" >&2
      exit 1
    fi
    analysis_files+=("../$anl_file")
  fi
done

if (( ${#analysis_files[@]} == 0 )); then
  echo "No supported GSI analysis diagnostics for $cdate in $PWD" >&2
  exit 1
fi

# Use the active Python environment when available. Otherwise use the RDASApp
# modules from the existing converter, with RDASApp available as an override.
if [[ -n ${RDASApp:-} ]] || ! python -c 'import netCDF4, numpy' >/dev/null 2>&1; then
  RDASApp=${RDASApp:-/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/PRs/rrfs-workflow.20260821.rw/sorc/RDASApp}
  if [[ ! -r "$RDASApp/ush/detect_machine.sh" ]]; then
    echo "Set RDASApp to your RDASApp build or activate Python with netCDF4 and numpy." >&2
    exit 1
  fi
  source "$RDASApp/ush/detect_machine.sh"
  module purge
  module use "$RDASApp/modulefiles"
  module load "RDAS/${MACHINE_ID}.intel"
fi

# Keep the converter's intermediate jdiag files isolated until it succeeds.
# Relative input paths keep "_anl" in an experiment name out of the lookup.
printf 'Converting GSI diagnostics for %s in %s\n' "$cdate" "$PWD"
work_dir=$(mktemp -d "$PWD/.convert_diags.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT
cd -- "$work_dir"
python "$script_dir/gdiag_to_jdiag.py" --surface-at2m "$cdate" "${analysis_files[@]}"

converted_files=(jdiag_*.nc)
if (( ${#converted_files[@]} == 0 )); then
  echo "The converter produced no diagnostics for $cdate." >&2
  exit 1
fi

for file in "${converted_files[@]}"; do
  output="gdiag_${file#jdiag_}"
  mv -f -- "$file" "../$output"
  # Replace any old temperature or humidity alias for this exact diagnostic.
  case "$output" in
    gdiag_*_airTemperatureAt2M_*|gdiag_*_specificHumidityAt2M_*)
      rm -f -- "../${output/At2M_/_}"
      ;;
    gdiag_*_airTemperature_*|gdiag_*_specificHumidity_*)
      stem=${output%_*}
      rm -f -- "../${stem}At2M_${output##*_}"
      ;;
  esac
  printf '  Wrote %s\n' "$output"
done

