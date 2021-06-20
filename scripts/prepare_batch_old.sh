#!/usr/bin/env bash
# Usage example:
#   prepare_batch.sh /path/to/instance/list.instances /path/to/results/dir [test / preprocess / bb / bb0]

if [ -z "$PROJ_DIR" ]
then
  if [ ! -z "${REPOS_DIR}" ]
  then
    echo "Please define PROJ_DIR (the root project dir, possibly ${REPOS_DIR}/vpc):"
  else
    echo "Please define PROJ_DIR (the root project dir):"
  fi
  read PROJ_DIR
  if [ -z "$PROJ_DIR" ]
    then echo "Need to define PROJ_DIR. Exiting."
    exit 1
  fi
fi

SILENT=1
MODE=preprocess
#MODE=bb
#MODE=bb0

export PROJ_DIR=`realpath -s ${PROJ_DIR}`
export VPC_DIR=${PROJ_DIR}

export INSTANCE_DIR=${PROJ_DIR}/data/instances
export RESULTS_DIR=${PROJ_DIR}/results
export RESULTS_DIR=/local1/$USER/results
export SCRIPT_DIR=${PROJ_DIR}/scripts
export SOL_DIR=${PROJ_DIR}/data/solutions

if [ $MODE == preprocess ]; then
  INSTANCE_LIST=${SCRIPT_DIR}/original.instances
else
  INSTANCE_LIST=${SCRIPT_DIR}/presolved.instances
fi

if [ ! -z $1 ]; then
  MODE=$1
fi
JOB_LIST="job_list_${MODE}.txt"

TASK_ID=0
> $JOB_LIST
while read line; do
  TASK_ID=$((TASK_ID+1))

  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi

  # Prepare out directory, based on current date
  CASE_NUM=`printf %03d $TASK_ID`
  STUB=`date +%F`
  OUT_DIR=${RESULTS_DIR}/$STUB/${MODE}/${CASE_NUM}

  # Print status (in silent mode, print a ".")
  if [ $SILENT != 1 ]; then
    echo "Preparing command to run instance $line (task $TASK_ID) at `date`"
  else
    echo -n "."
  fi

  # Check if solution exists
  arrIN=(${line//\// })
  arrIN[0]=""
  SOLFILE="$SOL_DIR"
  for entry in ${arrIN[@]}; do
    SOLFILE="${SOLFILE}/${entry}"
  done
  SOLFILE="${SOLFILE}_gurobi.mst.gz"
  if test -f "$SOLFILE"; then
    if [ $SILENT != 1 ]; then
      echo "$SOLFILE exists"
    fi
    SOLPARAM="--solfile=${SOLFILE}"
  else
    echo "*** WARNING: Could not find $SOLFILE"
    SOLPARAM=""
  fi

  # Finally, write command we will call to a file
  echo -n "mkdir -p ${OUT_DIR}; " >> ${JOB_LIST}
  if [ ${MODE} == preprocess ]; then
    echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps ${OUT_DIR} ${MODE} \" --temp=32 ${SOLPARAM}\" >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
  else
    echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps ${OUT_DIR} ${MODE} \" ${SOLPARAM}\" >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
  fi
done < ${INSTANCE_LIST}

echo "Done preparing $JOB_LIST"
