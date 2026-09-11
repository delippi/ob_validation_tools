#!/bin/bash
set -x
# This script is the main driver for converting GSI-diagnostic (gdiag) files
# to JEDI-diagnostic (jdiag) files. It processes gdiags from a single cycle
# and formats (both data and directory structure) to match JEDI to be used
# with the JEDI DA monitoring tools.
#
# Supported functions:
# - CONVERT_GDIAG_TO_JDIAG  : Converts gdiags to jdiags and saves the output in RRFSv2
#                   directory structure
#
#### USER-DEFINED VARIABLES #################################################
# Specify which functions to run (uncomment/comment to turn on/off)
CONVERT_GDIAG_TO_JDIAG="YES"

# Cycle start and end dates to process
SDATE=2024050601
EDATE=2024050601

# V2 retro experiment details (similar to rrfs-workflow/workflow/exp.setup)
VERSION="v1.2.0"
EXP_NAME="det_c13_3dvar_jedi_test"
if [[ $# -ge 2 ]]; then
  EXP_NAME="$1"
  SDATE="$2"
  EDATE=$SDATE
fi
sleep 5


OPSROOT="/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/MATS/ob_validation_tools/hofx_matcher/May2024_retro/${VERSION}"
COMROOT="${OPSROOT}/${EXP_NAME}"
JDIAGDIR="${COMROOT}"

# RRFSv1 EXP_NAME (RRFSv1 = benchmark)
CTL_VERSION="v1.2.0"
CTL_NAME="$EXP_NAME"
CTL_OPSROOT="${OPSROOT}"
CTL_COMROOT="${OPSROOT}/${CTL_NAME}"
CTL_JDIAGDIR="${CTL_COMROOT}"

# Specify your RDASApp build (mostly for module loads)
AppPath="/lfs/h2/emc/da/noscrub/donald.e.lippi/RRFSv1/PRs/rrfs-workflow.20260821.rw/sorc/"
RDASApp="RDASApp"

RDASApp="$AppPath/$RDASApp"
#### END OF USER-DEFINED VARIABLES ##########################################

# Start main execution
START=$(date +%s)

# Stop if the RDASApp directory does not exist.
if [[ ! -d "$RDASApp" ]]; then
  echo "Error: RDASApp directory does not exist: $RDASApp" >&2
  exit 1
fi

# Stop if the OPSROOT directory does not exist.
if [[ ! -d "$OPSROOT" ]]; then
  echo "Error: OPSROOT directory does not exist: $OPSROOT" >&2
  exit 1
fi

# Detect machine
source ${RDASApp}/ush/detect_machine.sh || exit 2

# Load necessary environment
module purge
module use ${RDASApp}/modulefiles || exit 3
module load RDAS/${MACHINE_ID}.intel || exit 4
export ndate=$(which ndate) || exit 5 # Load ndate
#module purge
#module load EVA/${MACHINE_ID}




if [[ -z "$ndate" ]]; then
  echo "Error: ndate command not found. Please ensure it is installed and available in your PATH." >&2
  exit 1
fi

# Initialize date with start date
date=${SDATE}

# Loop over dates from start to and including end date
while [[ ${date} -le ${EDATE} ]]; do
  pdy=${date:0:8}
  cyc=${date:8:10}

  if [[ ${CONVERT_GDIAG_TO_JDIAG:=NO} == "YES" ]] ; then
    echo "? Working on convert gdiag to jdiag: ${pdy} ${cyc}Z"
    # Array of full file paths (only need to specify _anl files)
    jdiags=(${JDIAGDIR}/rrfs.${pdy}/${cyc}/diag_conv_t_anl.${pdy}${cyc}.nc4)
    jdiags+=(${JDIAGDIR}/rrfs.${pdy}/${cyc}/diag_conv_q_anl.${pdy}${cyc}.nc4)
    jdiags+=(${JDIAGDIR}/rrfs.${pdy}/${cyc}/diag_conv_ps_anl.${pdy}${cyc}.nc4)
    jdiags+=(${JDIAGDIR}/rrfs.${pdy}/${cyc}/diag_conv_uv_anl.${pdy}${cyc}.nc4)
    jdiags+=(${JDIAGDIR}/rrfs.${pdy}/${cyc}/diag_conv_rw_anl.${pdy}${cyc}.nc4)

    # Run the Python script with analysis time and file list
    python gdiag_to_jdiag.py "$date" "${jdiags[@]}"

    for f in jdiag_*.nc; do
      [ -f "$f" ] || continue
      new="${f/jdiag_/gdiag_}"
      echo "Renaming $f -> $new"
      mv "$f" "$new"
    done

    mv gdiag*.nc ${JDIAGDIR}/rrfs.${pdy}/${cyc}/.

  fi

  # Increase date by 1 hour
  date=$(${ndate} 1 ${date})
done # date loop

exit 0
