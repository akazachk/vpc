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

export SCRIPT_DIR=${PROJ_DIR}/scripts
export INSTANCE_DIR=${PROJ_DIR}/data/instances
export INSTANCE_LIST=${SCRIPT_DIR}/slurm/small_original.instances
export RESULTS_DIR=${PROJ_DIR}/results

TASK_ID=0
while read line; do
  TASK_ID=$((TASK_ID+1))

  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi

  echo "Preparing command to run instance $line (task $TASK_ID) at `date`"
  echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps $RESULTS_DIR/bb/${TASK_ID} bb 2>&1" >> job_list_bb.txt
  echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps $RESULTS_DIR/bb0/${TASK_ID} bb0 2>&1" >> job_list_bb0.txt
  echo "nohup /usr/bin/time -v ${SCRIPT_DIR}/run_experiments.sh ${INSTANCE_DIR}/$line.mps $RESULTS_DIR/preprocess/${TASK_ID} preprocess 2>&1" >> job_list_preprocess.txt
done < ${INSTANCE_LIST}

