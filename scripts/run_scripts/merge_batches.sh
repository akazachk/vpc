#!/bin/bash
# The first argument is the results directory
# The second argument takes value for indicating the filename contained in each batch: "best", "test", "bb", "bb0", or "preprocess"
# This means we are merging vpc-best.csv, vpc-test.csv, or vpc-bb.csv from the batches, or cleaning_log.csv

MASTER_RESULTS_DIR="${VPC_DIR}/results"
SCRIPT_DIR="${VPC_DIR}/scripts/run_scripts"
EXECUTABLE="merge_results.sh"

if [ -z "$1" ]
  then RESULTS_DIR="${MASTER_RESULTS_DIR}/batches/test"
else
  RESULTS_DIR="$1"
fi

TMPNAME="vpc-bb" # default
if [ ! -z "$2" ]
then
  TMPNAME="vpc-${2}.csv"
  if [ "$2" = "preprocess" ]
    then TMPNAME="cleaning_log.csv"
  fi
fi

OUT_DIR="${RESULTS_DIR}"
OUTNAME="${OUT_DIR}/${TMPNAME}"

i=0
for batchname in `ls -d ${RESULTS_DIR}/*/`; do
  i=$(( $i + 1 ))  # maintain line count
  if [ "$2" = "full" ]
    then ${SCRIPT_DIR}/${EXECUTABLE} ${batchname} 0 full
  fi
  if [ $i = 1 ]
    then 
      head -n 2 ${batchname}${TMPNAME} > ${OUTNAME}
  fi
  echo "Copying ${TMPNAME} from ${batchname}"
  #tail -n +3 ${batchname}${TMPNAME} >> ${OUTNAME}
  tail -n +3 ${batchname}${TMPNAME} | grep DONE >> ${OUTNAME}
done

echo "Copying and sorting merged batch result from ${OUTNAME} to ${MASTER_RESULTS_DIR}/${TMPNAME}"
(head -n 2 ${OUTNAME} && (tail -n +3 ${OUTNAME} | sort)) > ${MASTER_RESULTS_DIR}/${TMPNAME}
cp ${MASTER_RESULTS_DIR}/${TMPNAME} ${OUTNAME}
