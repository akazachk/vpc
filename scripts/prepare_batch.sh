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
export OPTFILE="${VPC_DIR}/data/ip_obj.csv"
export RESULTS_DIR=${PROJ_DIR}/results
export SCRIPT_DIR=${PROJ_DIR}/scripts
export SOL_DIR=${PROJ_DIR}/data/solutions

export INSTANCE_DIR=/local1/$USER/instances
export RESULTS_DIR=/local1/$USER/results

EXECUTABLE="${PROJ_DIR}/Release/vpc"

if [ $MODE == preprocess ]; then
  INSTANCE_LIST=${SCRIPT_DIR}/original.instances
else
  INSTANCE_LIST=${SCRIPT_DIR}/presolved.instances
fi

# Accept user options for instance list, results directory, and mode
if [ ! -z $1 ]; then
  INSTANCE_LIST=$1
fi
if [ ! -z $2 ]; then
  RESULTS_DIR=$2
fi
if [ ! -z $3 ]; then
  MODE=$3
fi
JOB_LIST="job_list_${MODE}.txt"

# Set parameters
PARAMS=" --optfile=${OPTFILE}"
if [ $MODE == bb ]; then
  depthList=(2 4 8 16 32 64)
  PARAMS="$PARAMS -t 3600"
  PARAMS="$PARAMS --rounds=1"
  PARAMS="$PARAMS --bb_runs=1"
  PARAMS="$PARAMS --bb_mode=10"
  PARAMS="$PARAMS --bb_timelimit=3600"
  PARAMS="$PARAMS --use_all_ones=1"
  PARAMS="$PARAMS --use_iter_bilinear=1"
  PARAMS="$PARAMS --use_disj_lb=1"
  PARAMS="$PARAMS --use_tight_points=0"
  PARAMS="$PARAMS --use_tight_rays=0"
  PARAMS="$PARAMS --use_unit_vectors=0"
  PARAMS="$PARAMS --gomory=-1"
elif [ $MODE == bb0 ]; then
  depthList=(0)
  PARAMS="$PARAMS -t 3600"
  PARAMS="$PARAMS --rounds=1"
  PARAMS="$PARAMS --bb_runs=7"
  PARAMS="$PARAMS --bb_mode=001"
  PARAMS="$PARAMS --bb_timelimit=3600"
elif [ $MODE == preprocess ]; then
  depthList=(0)
  PARAMS="$PARAMS -t 7200"
  PARAMS="$PARAMS --preprocess=1"
  PARAMS="$PARAMS --bb_runs=1"
  PARAMS="$PARAMS --bb_mode=001"
  # Set strategy (536 uses Gurobi, 532 uses CPLEX)
  PARAMS="$PARAMS --bb_strategy=536"
  PARAMS="$PARAMS --bb_timelimit=7200"
  PARAMS="$PARAMS --temp=32"
elif [ $MODE == test ]; then
  depthList=(2)
else
  echo "*** ERROR: Option $MODE not recognized"
  exit
fi

TASK_ID=0
TOTAL_ERRORS=0
> $JOB_LIST
for d in ${depthList[*]}; do
  echo "Depth $d"
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

    # Check if file exists
    FILE=${INSTANCE_DIR}/$line.mps
    if [ ! -f "$FILE" ] && [ ! -f "$FILE.gz" ] && [ ! -f "$FILE.bz2" ]; then
      echo "*** ERROR: $FILE does not exist; skipping."
      TOTAL_ERRORS=$((TOTAL_ERRORS+1))
    else
      # Finally, write command we will call to a file
      echo -n "mkdir -p ${OUT_DIR}; " >> ${JOB_LIST}
      echo "nohup /usr/bin/time -v $EXECUTABLE -f ${FILE} --logfile=${OUT_DIR}/vpc-${MODE}.csv $SOLPARAM $PARAMS -d$d >> ${OUT_DIR}/log.out 2>&1" >> ${JOB_LIST}
    fi
  done < ${INSTANCE_LIST}
done # loop over depth list

echo "Done preparing $JOB_LIST"
