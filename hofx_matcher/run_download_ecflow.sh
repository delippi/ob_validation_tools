#!/bin/bash

set -euo pipefail
shopt -s nullglob

usage() {
  printf 'Usage: %s EXP_NAME YYYYMMDDHH [JOB_ID.BATCH_HOST]\n' "$0" >&2
  printf '       %s --pair GSI_EXP_NAME JEDI_EXP_NAME YYYYMMDDHH\n' "$0" >&2
  exit 2
}

# ecFlow runs GSI and JEDI in separate experiments. Reuse the single-run
# download for each experiment so each has its own menu and staging directory.
if [[ ${1:-} == --pair ]]; then
  [[ $# -eq 4 ]] || usage
  if [[ ( $2 != *_gsi_* && $2 != *_gsi ) ||
        ( $3 != *_jedi_* && $3 != *_jedi ) ]]; then
    echo "Pass the GSI experiment first and the JEDI experiment second." >&2
    exit 2
  fi
  "$BASH" "${BASH_SOURCE[0]}" "$2" "$4"
  "$BASH" "${BASH_SOURCE[0]}" "$3" "$4"
  exit 0
fi

# Run from hofx_matcher. Omit the job suffix to choose from available runs.
[[ $# -ge 2 && $# -le 3 ]] || usage

expname=$1
cdate=$2
run_id_pattern='^[0-9]+\.[[:alnum:]_][[:alnum:]_.-]*$'

if [[ ! $expname =~ ^[[:alnum:]_][[:alnum:]_.-]*$ ||
      ! $cdate =~ ^[0-9]{8}([01][0-9]|2[0-3])$ ]]; then
  echo "Invalid experiment name or cycle." >&2
  exit 2
fi

case "$expname" in
  *_jedi_*|*_jedi) analysis=jedi ;;
  *_gsi_*|*_gsi) analysis=gsi ;;
  *)
    echo "Expected a jedi or gsi component in the experiment name: $expname" >&2
    exit 2
    ;;
esac

pdy=${cdate:0:8}
cyc=${cdate:8:2}
source_root=${ECFLOW_ROOT:-/lfs/h2/emc/stmp/${USER}/ecflow_rrfs}
experiment_dir="$source_root/$expname/$pdy/$expname"

# GSI diagnostics are produced by a separate task with its own job ID.
task="analysis_${analysis}"
if [[ $analysis == gsi ]]; then
  task=analysis_gsi_diag
fi
run_prefix="rrfs_det_${task}_${cyc}"
analysis_dir="$experiment_dir/rrfs_${task}_${cyc}_v1.0_prod/det"

# Some layouts keep the diagnostic job directory under the GSI analysis task.
# Keep the diagnostic job prefix in either location.
if [[ $analysis == gsi && ! -d $analysis_dir ]]; then
  gsi_analysis_dir="$experiment_dir/rrfs_analysis_gsi_${cyc}_v1.0_prod/det"
  if [[ -d $gsi_analysis_dir ]]; then
    analysis_dir=$gsi_analysis_dir
  fi
fi

if [[ $# -eq 3 ]]; then
  run_id=$3
  if [[ ! $run_id =~ $run_id_pattern ]]; then
    echo "Invalid job suffix: $run_id" >&2
    exit 2
  fi
else
  # List only run directories for this analysis and cycle.
  run_ids=()
  for run_dir in "$analysis_dir/$run_prefix".*; do
    suffix=${run_dir##*/$run_prefix.}
    if [[ -d $run_dir && $suffix =~ $run_id_pattern ]]; then
      run_ids+=("$suffix")
    fi
  done

  if (( ${#run_ids[@]} == 0 )); then
    echo "No analysis run directories found in: $analysis_dir" >&2
    exit 1
  fi

  # Sort on the numeric job ID so the earliest submission is option 1.
  mapfile -t run_ids < <(printf '%s\n' "${run_ids[@]}" | LC_ALL=C sort -t . -k1,1n)
  printf '%s jobs for %s %s (ascending job ID)\n' "$task" "$expname" "$cdate" >&2
  PS3='Select a run number: '
  run_id=
  select run_id in "${run_ids[@]}"; do
    if [[ -n $run_id ]]; then
      break
    fi
    echo "Choose a run number from the list." >&2
  done

  if [[ -z $run_id ]]; then
    echo "No run selected." >&2
    exit 1
  fi
fi

run_name="$run_prefix.$run_id"
source_dir="$analysis_dir/$run_name"
dest_dir="ecflow/$expname/rrfs.$pdy/$cyc/$run_name"

if [[ ! -d $source_dir ]]; then
  echo "Analysis run directory does not exist: $source_dir" >&2
  exit 1
fi

# Copy only the diagnostics used by the matcher and converter. GSI can
# leave compressed NetCDF files in the work directory after concatenation.
if [[ $analysis == jedi ]]; then
  diagnostics=("$source_dir"/jdiag_*.nc)
else
  diagnostics=("$source_dir"/diag_*."$cdate".nc4 "$source_dir"/diag_*."$cdate".nc4.gz)
fi

if (( ${#diagnostics[@]} == 0 )); then
  echo "No $analysis NetCDF diagnostics found in: $source_dir" >&2
  exit 1
fi

mkdir -p "$dest_dir"
printf 'Staging %s\n     to %s\n' "$source_dir" "$dest_dir"
rsync -av -- "${diagnostics[@]}" "$dest_dir/"

# Decompress the staged copies so gdiag_to_jdiag.py can open them directly.
for file in "${diagnostics[@]}"; do
  if [[ $file == *.nc4.gz ]]; then
    gunzip -f -- "$dest_dir/${file##*/}"
  fi
done

