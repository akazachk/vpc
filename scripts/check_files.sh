#!/usr/bin/env bash
# Usage:
#   check_files.sh presolved.instances /path/to/instances/
#
# Argument 1: instance list
# Argument 2: instance directory

check_proj_dir() {
  if [ -z "$PROJ_DIR" ]
  then
    if [ ! -z "${REPOS_DIR}" ]
    then
      echo "Please define PROJ_DIR (the root project dir, possibly ${REPOS_DIR}/vpc):" > /dev/stderr
    else
      echo "Please define PROJ_DIR (the root project dir):" > /dev/stderr
    fi
    read PROJ_DIR
  fi
  echo $PROJ_DIR
}

if [ ! -z $1 ]; then
  INSTANCE_LIST="$1"
else
  PROJ_DIR=$(check_proj_dir)
  if [ -z "$PROJ_DIR" ]
    then echo "Need to define PROJ_DIR. Exiting." > /dev/stderr
    exit 1
  fi
  SCRIPT_DIR=${PROJ_DIR}/scripts
  INSTANCE_LIST=${SCRIPT_DIR}/small_presolved.instances
fi

if [ ! -z $2 ]; then
  INSTANCE_DIR="$2"
else
  PROJ_DIR=$(check_proj_dir)
  if [ -z "$PROJ_DIR" ]
    then echo "Need to define PROJ_DIR. Exiting." > /dev/stderr
    exit 1
  fi
  INSTANCE_DIR=${PROJ_DIR}/data/instances
fi

echo "Instance directory set to $INSTANCE_DIR"
echo "Instance list is $INSTANCE_LIST"

FOUND=0
TOTAL=0
MISSING="Missing:"
while read line; do
  TOTAL=$((TOTAL+1))
  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi

  FILE=${INSTANCE_DIR}/$line.mps
  if [ ! -f "$FILE" ] && [ ! -f "$FILE.gz" ] && [ ! -f "$FILE.bz2" ]; then
    echo "$FILE does not exist"
    MISSING="$MISSING\n$line"
  else
    echo -n "."
    FOUND=$((FOUND+1))
  fi
done < ${INSTANCE_LIST}

echo ""
echo "Done! Found $FOUND/$TOTAL files."
echo ""
echo -e $MISSING
