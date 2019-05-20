#!/bin/bash
# VPC_DIR should be defined for these scripts
#
# Three arguments: 
# 1. 'bb' or another suffix (sets script to be used as vpc/run_vpc_$1.py) [optional, but if this is not used, then none of subsequent options can be given; default: bb]
# 2. full path to instance list, which should be located in the same directory as the instances [optional; default: ${VPC_DIR}/data/instances/test.instances]
# 3. '1' or a '0', referring to whether batch mode is used [optional; default: 0]

# Defaults
export VPC_DIR=${HOME}/repos/vpc
DEFAULT_INSTANCE_DIR="${VPC_DIR}/data/instances"
DEFAULT_INSTANCE_LIST="${VPC_DIR}/data/instances/test.instances"
DEFAULT_BATCH_LIST="${VPC_DIR}/data/instances/test.batch"
INSTANCE_LIST="${DEFAULT_INSTANCE_LIST}"
BATCH_MODE=0

SCRIPT_DIR="${VPC_DIR}/scripts/run_scripts"
OUT_DIR="${VPC_DIR}/results"

DEBUG="echo"
DEBUG=""

CUT_TYPE="vpc"

# Process run type
if [ -z "$1" ]
then 
  #echo "*** ERROR: Need to specify run type."
  #exit 1
  export RUN_TYPE_STUB="bb"
else
  export RUN_TYPE_STUB="$1"
fi

# Process other arguments
if [ -z "$2" ]
then 
  INSTANCE_LIST="${DEFAULT_INSTANCE_LIST}"
elif [ "$2" = 1 ]
then 
  BATCH_MODE=1
  INSTANCE_LIST="${DEFAULT_BATCH_LIST}"
elif [ "$2" = 0 ]
then 
  BATCH_MODE=0
  INSTANCE_LIST="${DEFAULT_BATCH_LIST}"
else 
  INSTANCE_LIST="$2"
  if [ "$3" = 1 ]
    then 
      BATCH_MODE=1
  elif [ "$3" = 0 ]
    then 
      BATCH_MODE=0
  fi
fi

TMPNAME="run_${CUT_TYPE}_${RUN_TYPE_STUB}"
TMPNAME_EXT=".py"

echo "Options selected: run ${CUT_TYPE}/${TMPNAME} using instances given in ${INSTANCE_LIST}. In batch mode? ${BATCH_MODE}."
#exit 1

# Proceed depending on whether run is in batches or not
if [ $BATCH_MODE == 0 ] 
then
  echo "Running ${CUT_TYPE}/${TMPNAME}${TMPNAME_EXT} from ${INSTANCE_LIST}"
  nohup python ${SCRIPT_DIR}/${CUT_TYPE}/${TMPNAME}${TMPNAME_EXT} 0 ${INSTANCE_LIST} >& ${OUT_DIR}/nohup.out &
else
  FSTUB=`head -n 1 ${INSTANCE_LIST}`
  tmpfilename=""
  for line in `tail -n +2 ${INSTANCE_LIST}`; do
    # Skip empty lines
    if [ -z "$line" ]
    then
      continue
    fi
    
    #if [ ! -z $line ]
    #then
    #  echo "Current line: $line"
    #  echo "Current batch: $batch"
    #fi

    # If line ends with '/', then we start a new batch
    len="$((${#line}-1))"
    if [ "${line:len:1}" == "/" ]
    then
      # If there was an old batch, then run it
      if [ ! -z "${tmpfilename}" ]
      then
        echo "Running ${CUT_TYPE}/${TMPNAME}${TMPNAME_EXT} from ${tmpfilename}"
        nohup python ${SCRIPT_DIR}/${CUT_TYPE}/${TMPNAME}${TMPNAME_EXT} 1 ${tmpfilename} >& ${OUT_DIR}/nohup.out &
      fi  

      # Now we create the new batch
      batchstub="${line:0:len}"
      tmpfilename="/tmp/${FSTUB}.batch${batchstub}XXX"
      tmpfilename=$(mktemp -q ${tmpfilename})
      echo "${FSTUB}" > "${tmpfilename}"
      echo "${line}" >> "${tmpfilename}"
    else
      # Else, we add the current line to the current batch
      echo "${line}" >> "${tmpfilename}"
    fi  
  done

  # Now process the last batch
  if [ ! -z "${tmpfilename}" ]
  then
    echo "Running ${CUT_TYPE}/${TMPNAME}${TMPNAME_EXT} from ${tmpfilename}"
    nohup python ${SCRIPT_DIR}/${CUT_TYPE}/${TMPNAME}${TMPNAME_EXT} 1 ${tmpfilename} >& ${OUT_DIR}/nohup.out &
  fi  
fi
