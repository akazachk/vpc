#!/usr/bin/env bash
# Example of calling this script:
# ./run_experiments.sh test.instances ../results [test, bb, bb0, preprocess,...]
#
# The first argument is either a list of instances (each instance is assumed to be located in ${INSTANCE_DIR}, defined below as ${VPC_DIR}/data/instances) or the full path to an instance (lp/mps file)
# The second argument is where results will go
# The third argument is 'bb' or another suffix (sets script to be used as python/run_vpc_$1.py) [optional, but if this is not used, then subsequent options cannot be given; default: bb]
#
# If a list of instances is given, it must either be with extension ".instances" or ".batch"
# The former indicates that the instances should be run sequentially
# A ".batch" file will have batches of instances, where the instances within each batch will be run sequentially
# The first line of this file will be the directory "stub"
# and output will be sent to ${VPC_DIR}/results/stub if batch mode is off,
# and to ${VPC_DIR}/results/batches/stub/batchname if it is on.
# Each line contains either a relative input path to an instance (with or without the extension) or a batch name.
# The path is relative to ${VPC_DIR}/data/instances, e.g., the line will be original/miplib2/bm23.
# A batch name is distinguished by having the batch end with a '/', e.g., '2/' or 'batch2/'.

# Defaults
if [ -z "$VPC_DIR" ]
then 
  if [ ! -z "${REPOS_DIR}" ]
  then
    echo "Please define VPC_DIR (the root vpc dir, possibly ${REPOS_DIR}/vpc):"
  else
    echo "Please define VPC_DIR (the root vpc dir):"
  fi
  read VPC_DIR
  if [ -z "$VPC_DIR" ]
    then echo "Need to define VPC_DIR. Exiting."
    exit 1
  fi
fi
export VPC_DIR=${VPC_DIR} # used in python script
export INSTANCE_DIR="${VPC_DIR}/data/instances" # used in python script
SCRIPT_DIR="${VPC_DIR}/scripts"
INSTANCE_LIST="${VPC_DIR}/scripts/test.batch"
BATCH_MODE=0
CUT_TYPE="vpc"

# Process instance list
if [ -z "$1" ]
then
  echo "*** ERROR: Need to provide instance list as first argument."
  exit 1
else
  INSTANCE_LIST="$1"
fi

# Process results directory
if [ -z "$2" ]
then
  echo "*** ERROR: Need to provide results directory as second argument."
  exit 1
else
  RESULTS_DIR="$2"
fi

# Process run type
if [ -z "$3" ]
then 
  #echo "*** ERROR: Need to specify run type."
  #exit 1
  export RUN_TYPE_STUB="bb"
else
  export RUN_TYPE_STUB="$3"
fi

#export INSTANCE_DIR=${INSTANCE_LIST%/*}

SCRIPTNAME="run_${CUT_TYPE}.py"

echo "Running experiments from ${SCRIPT_DIR}/${SCRIPTNAME} with instances in list ${INSTANCE_LIST} and instance dir = ${INSTANCE_DIR}."

# Proceed depending on whether run is in batches or not
line=${INSTANCE_LIST}
tmplenlp="$((${#line}-3))"
tmplenmps="$((${#line}-4))"
tmplenlpgz="$((${#line}-6))"
tmplenmpsgz="$((${#line}-7))"
tmplenlpbz="$((${#line}-7))"
tmplenmpsbz="$((${#line}-8))"
tmplentxt="$((${#line}-4))"
tmplenbatch="$((${#line}-6))"
tmpleninst="$((${#line}-10))"
if [ "${line:$tmpleninst:10}" == ".instances" ] || [ "${line:$tmplenlp:3}" == ".lp" ] || [ "${line:$tmplenmps:4}" == ".mps" ] || [ "${line:$tmplenlpgz:6}" == ".lp.gz" ] || [ "${line:$tmplenmpsgz:7}" == ".mps.gz" ] || [ "${line:$tmplenlpbz:7}" == ".lp.bz2" ] || [ "${line:$tmplenmpsbz:8}" == ".mps.bz2" ]
then
  echo "Running ${SCRIPT_DIR}/${SCRIPTNAME} from ${INSTANCE_LIST} in sequential mode, output sent to ${RESULTS_DIR}"
  nohup python ${SCRIPT_DIR}/${SCRIPTNAME} ${RUN_TYPE_STUB} ${INSTANCE_LIST} ${RESULTS_DIR} >& ${RESULTS_DIR}/nohup.out &
elif [ "${line:$tmplenbatch:6}" == ".batch" ]
then
  echo "Running ${SCRIPT_DIR}/${SCRIPTNAME} from ${INSTANCE_LIST} in batch mode, output sent to ${RESULTS_DIR}"
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
        echo "Running ${SCRIPT_DIR}/${SCRIPTNAME} from ${tmpfilename}"
        nohup python ${SCRIPT_DIR}/${SCRIPTNAME} ${RUN_TYPE_STUB} ${tmpfilename} ${RESULTS_DIR} >& ${RESULTS_DIR}/nohup.out &
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
    echo "Running ${SCRIPT_DIR}/${SCRIPTNAME} from ${tmpfilename}"
    nohup python ${SCRIPT_DIR}/${SCRIPTNAME} ${RUN_TYPE_STUB} ${tmpfilename} ${RESULTS_DIR} >& ${RESULTS_DIR}/nohup.out &
  fi  
else
  echo "Could not identify type of instance file given by $line"
fi
