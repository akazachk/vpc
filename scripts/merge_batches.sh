#!/usr/bin/env bash
# The first argument is the results directory
# The second argument takes value for indicating the filename contained in each batch: "best", "test", "bb", "bb0", or "preprocess"
# This means we are merging vpc-test.csv, vpc-bb.csv, etc. from the batches, or cleaning_log.csv

if [ -z "$1" ]
then
  echo "Need to specify where the results are as the first argument. Exiting."
  exit
else
  RESULTS_DIR="$1"
fi

TMPNAME="vpc-bb" # default
if [ -z "$2" ]
then
  echo "Need to define the filename stub as the second argument. Exiting."
  exit
else
  TMPNAME="vpc-${2}.csv"
  if [ "$2" = "preprocess" ]
    then TMPNAME="cleaning_log.csv"
  fi
fi

OUTNAME="${RESULTS_DIR}/${TMPNAME}"
ERR_OUTNAME="${RESULTS_DIR}/errors_${TMPNAME}"

i=0
for batchname in `ls -d ${RESULTS_DIR}/*/`; do
  i=$(( $i + 1 ))  # maintain line count
  if [ $i = 1 ]
    then
      head -n 2 ${batchname}${TMPNAME} > ${OUTNAME}
      head -n 2 ${batchname}${TMPNAME} > ${ERR_OUTNAME}
  fi

  echo "Copying ${TMPNAME} from ${batchname} to ${OUTNAME}"
  #tail -n +3 ${batchname}${TMPNAME} >> ${OUTNAME}
  tail -n +3 ${batchname}${TMPNAME} | grep DONE >> ${OUTNAME}

  echo "Copying errors in ${TMPNAME} from ${batchname} to ${ERR_OUTNAME}"
  tail -n +3 ${batchname}${TMPNAME} | grep "ERROR" >> ${ERR_OUTNAME}
done

#if [ ! -z $BASH_VERSION ]
#then
#    SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
#elif [ ! -z $ZSH_VERSION ]
#then
#    SCRIPT_DIR="${0:a:h}"
#else
#    SCRIPT_DIR="${0%/*}"
#fi
#MASTER_RESULTS_DIR="${SCRIPT_DIR}/results_test"
#echo "Copying and sorting merged batch result from ${OUTNAME} to ${MASTER_RESULTS_DIR}/${TMPNAME}"
#(head -n 2 ${OUTNAME} && (tail -n +3 ${OUTNAME} | sort)) > ${MASTER_RESULTS_DIR}/${TMPNAME}
#cp ${MASTER_RESULTS_DIR}/${TMPNAME} ${OUTNAME}
