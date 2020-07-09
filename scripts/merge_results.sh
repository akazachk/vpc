#!/usr/bin/env bash
# Arguments: results directory, batch mode flag (0 or 1), run type (test, bb, bb0, or preprocess)
# Default is ${VPC_DIR}/results/test 0 test
#
# The first argument can be the (full) directory containing the results or batch folders containing results
# Alternatively, it can be a 1 (batch mode) or 0 (sequential mode)
# If the first argument is not 1, then the second argument may be a 1 indicating batch mode
# The argument after the 1 means the results are contained in vpc-test.csv, vpc-bb.csv, etc.

# Get the directory in which scripts are located (note that this will *fail* in some cases, such as when this script is sourced)
# We assume that the merge_batches script is located in the same directory as this one
if [ ! -z $BASH_VERSION ]
then
    SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
elif [ ! -z $ZSH_VERSION ]
then
    SCRIPT_DIR="${0:a:h}"
else
    SCRIPT_DIR="${0%/*}"
fi
MASTER_RESULTS_DIR="${SCRIPT_DIR}/results_test"
RESULTS_DIR="${MASTER_RESULTS_DIR}/test"
RUN_TYPE_STUB="test"
export CUT_TYPE="vpc"

# Set the RESULTS_DIR and run in batch mode if so specified
if [ "$1" = 1 -o "$2" = 1 ]
then
  # Batch mode
  if [ "$1" != "1" ]
  then
    # Second argument is a 1
    RESULTS_DIR="$1"
    RUN_TYPE_STUB="$3"
  else
    # First argument is a 1
    RUN_TYPE_STUB="$2"
  fi
  echo "Combining results from batches in ${RESULTS_DIR}"
  ${SCRIPT_DIR}/merge_batches.sh "${RESULTS_DIR}" "${RUN_TYPE_STUB}"
  exit 1
elif [ ! -z "$1" ]
then
  if [ "$1" != 0 -a "$1" != "best" -a "$1" != "test" -a "$1" != "bb" -a "$1" != "bb0" -a "$1" != "preprocess"]
  then
    RESULTS_DIR="$1"
  fi
fi

if [ "$1" = "test" -o "$2" = "test" -o "$3" = "test" ]
then
  RUN_TYPE_STUB="test"
elif [ "$1" = "best" -o "$2" = "best" -o "$3" = "best" ]
then
  RUN_TYPE_STUB="best"
elif [ "$1" = "bb" -o "$2" = "bb" -o "$3" = "bb" ]
then
  RUN_TYPE_STUB="bb"
elif [ "$1" = "bb0" -o "$2" = "bb0" -o "$3" = "bb0" ]
then
  RUN_TYPE_STUB="bb0"
elif [ "$1" = "preprocess" -o "$2" = "preprocess" -o "$3" = "preprocess" ]
then
  RUN_TYPE_STUB="preprocess"
fi

TMPNAME="${CUT_TYPE}-${RUN_TYPE_STUB}.csv"
if [ "$RUN_TYPE_STUB" = "preprocess" ]
  then TMPNAME="cleaning_log.csv"
fi

echo "Currently sorting results turned off"
#echo "Copying and sorting results from ${RESULTS_DIR}/${TMPNAME} to ${MASTER_RESULTS_DIR}/${TMPNAME}"
#(head -n 2 ${RESULTS_DIR}/${TMPNAME} && (tail -n +3 ${RESULTS_DIR}/${TMPNAME} | sort)) > ${MASTER_RESULTS_DIR}/${TMPNAME}
