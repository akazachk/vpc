#!/bin/bash

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

export PROJ_DIR=`realpath -s ${PROJ_DIR}`
export SCRIPT_DIR=${PROJ_DIR}/scripts
export INSTANCE_DIR=${PROJ_DIR}/data/instances
export INSTANCE_LIST=${SCRIPT_DIR}/small_presolved.instances
export RESULTS_DIR=${PROJ_DIR}/results

TASK_ID=0
> job_list_bb.txt
> job_list_bb0.txt
while read line; do
  TASK_ID=$((TASK_ID+1))

  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi

  CASE_NUM=`printf %03d $TASK_ID`
  STUB=`date +%F`
  echo "Preparing command to run instance $line (task $TASK_ID) at `date`"
  echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps $RESULTS_DIR/$STUB/bb/${CASE_NUM} bb 2>&1" >> job_list_bb.txt
  echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps $RESULTS_DIR/$STUB/bb0/${CASE_NUM} bb0 2>&1" >> job_list_bb0.txt
  #echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps $RESULTS_DIR/$STUB/preprocess/${CASE_NUM} preprocess 2>&1" >> job_list_preprocess.txt
done < ${INSTANCE_LIST}

