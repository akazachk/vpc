#!/usr/bin/env bash
# Arguments: results directory, run type (test, bb, bb0, or preprocess), [optional: batch mode flag (0 or 1)]
# Example: ./merge.sh ${VPC_DIR}/scripts/results_test test 1
#
# The first argument is the (full) directory containing the results or batch folders containing results
# The second argument indicates the results are contained in vpc-test.csv, vpc-bb.csv, etc.
# The last (optional) argument is a 1 (batch mode) or 0 (sequential mode); default is 1

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
RESULTS_DIR="${SCRIPT_DIR}/results_test"
BATCH_MODE=1
RUN_TYPE_STUB="test"
export CUT_TYPE="vpc"

# Set the RESULTS_DIR
if [ -z $1 ]; then
  echo "First argument must be results directory."
  exit 1
else
  RESULTS_DIR="$1"
fi

# Set the RUN_TYPE_STUB
if [ -z $2 ]; then
  echo "Second argument must be the stub (results from vpc-{stub}.csv will be merged)."
  exit 1
else
  RUN_TYPE_STUB="$2"
fi

# Set batch mode
if [ "$3" = 1 ]
then
  BATCH_MODE=1
elif [ "$3" = 0 ]; then
  BATCH_MODE=0
elif [ ! -z $3 ]; then
  echo "Third argument must be a 0 (sequential mode), a 1 (batch mode), or empty."
  exit 1
fi

TMPNAME="${CUT_TYPE}-${RUN_TYPE_STUB}.csv"
if [ "$RUN_TYPE_STUB" = "preprocess" ]
  then TMPNAME="cleaning_log.csv"
fi
OUTNAME="${RESULTS_DIR}/${TMPNAME}"
ERR_OUTNAME="${RESULTS_DIR}/errors_${TMPNAME}"

if [ $BATCH_MODE = 1 ]; then
  echo "Combining results from batches in ${RESULTS_DIR}"
  i=0
  for batchname in `ls -d ${RESULTS_DIR}/*/`; do
    i=$(( $i + 1 ))  # maintain line count
    if [ $i = 1 ]
      then
        head -n 2 ${batchname}${TMPNAME} > ${OUTNAME}
        head -n 2 ${batchname}${TMPNAME} > ${ERR_OUTNAME}
    fi

    echo "Copying ${TMPNAME} from ${batchname} to ${OUTNAME}"
    tail -n +3 ${batchname}${TMPNAME} | grep DONE >> ${OUTNAME}

    echo "Copying errors in ${TMPNAME} from ${batchname} to ${ERR_OUTNAME}"
    tail -n +3 ${batchname}${TMPNAME} | grep "ERROR" >> ${ERR_OUTNAME}
  done
else
  echo "Currently sorting results in sequential mode turned off"
  #echo "Copying and sorting results from ${RESULTS_DIR}/${TMPNAME} to ${MASTER_RESULTS_DIR}/${TMPNAME}"
  #(head -n 2 ${RESULTS_DIR}/${TMPNAME} && (tail -n +3 ${RESULTS_DIR}/${TMPNAME} | sort)) > ${MASTER_RESULTS_DIR}/${TMPNAME}
fi
