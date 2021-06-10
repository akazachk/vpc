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
export INSTANCE_LIST=${SCRIPT_DIR}/small_presolved.instances
export INSTANCE_DIR=${PROJ_DIR}/data/instances

if [ ! -z $1 ]; then
  INSTANCE_LIST="$1"
fi
if [ ! -z $2 ]; then
  INSTANCE_DIR="$2"
fi

echo "Instance directory set to $INSTANCE_DIR"
echo "Instance list is $INSTANCE_LIST"

FOUND=0
TOTAL=0
while read line; do
  TOTAL=$((TOTAL+1))
  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi

  FILE=${INSTANCE_DIR}/$line.mps
  if [ ! -f "$FILE" ] && [ ! -f "$FILE.gz" ] && [ ! -f "$FILE.bz2" ]; then
    echo "$FILE does not exist"
  else
    echo -n "."
    FOUND=$((FOUND+1))
  fi
done < ${INSTANCE_LIST}

echo ""
echo "Done! Found $FOUND/$TOTAL files."
