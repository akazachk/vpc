#!/usr/bin/env bash

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
export INSTANCE_LIST=${SCRIPT_DIR}/slurm/small_presolved.instances

while read line; do
  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi

  FILE=${INSTANCE_DIR}/$line.mps.gz
  if [ ! -f "$FILE" ]; then
    echo "$FILE does not exist"
  fi
done < ${INSTANCE_LIST}
