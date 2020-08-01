#!/usr/bin/env bash
# Example of calling this script:
# ./run_experiments.sh test.instances ../results/test [test, bb, bb0, preprocess,...] "extra_arg more_extra_args"
#
# The first argument is either a list of instances (each instance is assumed to be located in ${INSTANCE_DIR}, defined below as ${VPC_DIR}/data/instances) or the full path to an instance (lp/mps file)
# The second argument is where results will go, saved into RESULTS_DIR
# The third argument is 'bb' or another suffix (sets script to be used as python3/run_vpc_$1.py) [optional, but if this is not used, then subsequent options cannot be given; default: bb]
# The fourth argument is anything extra we wish to pass to the solver, passed in quotation marks so it just counted as one argument
#
# If a list of instances is given, it must either be with extension ".instances" or ".batch"
# The former indicates that the instances should be run sequentially
# A ".batch" file will have batches of instances, where the instances within each batch will be run sequentially
# Output will be sent to ${RESULTS_DIR} if batch mode is off, and to ${RESULTS_DIR}/batchname if it is on.
# Each line contains either a relative input path to an instance (with or without the extension) or a batch name.
# The path is relative to ${VPC_DIR}/data/instances, e.g., the line for bm23 from miplib2 will be original/miplib2/bm23, or presolved/miplib2/bm23_presolved for the presolved version.
# A batch name is distinguished by having the line end with a '/', e.g., '2/' or 'batch2/'.

split_on_commas() {
  local IFS=,
  local WORD_LIST=($1)
  for word in "${WORD_LIST[@]}"; do
    echo "$word"
  done
}

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
export VPC_DIR=${VPC_DIR} # used in python3 script
export INSTANCE_DIR="${VPC_DIR}/data/instances" # used in python3 script
SCRIPT_DIR="${VPC_DIR}/scripts"
INSTANCE_LIST="${VPC_DIR}/scripts/test.batch"
SCRIPTNAME="run_vpc.py"
CURRSCRIPTNAME=${0##*/}

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
  mkdir -p ${RESULTS_DIR}
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

# Disabled: Define the instance directory from which relative paths are given in the instance list
#export INSTANCE_DIR=${INSTANCE_LIST%/*}

#echo "Running experiments from ${SCRIPT_DIR}/${SCRIPTNAME} with instance list ${INSTANCE_LIST}, assuming instances are in ${INSTANCE_DIR}, output sent to ${RESULTS_DIR}."

# Proceed depending on whether run is in batches or not
line=${INSTANCE_LIST}
tmplenlp="$((${#line}-3))"
tmplenmps="$((${#line}-4))"
tmplenlpgz="$((${#line}-6))"
tmplenmpsgz="$((${#line}-7))"
tmplenlpbz="$((${#line}-7))"
tmplenmpsbz="$((${#line}-8))"
tmplentxt="$((${#line}-4))"
tmpext=${INSTANCE_LIST##*.}
if [ "${line:$tmplenlp:3}" == ".lp" ] || [ "${line:$tmplenmps:4}" == ".mps" ] || [ "${line:$tmplenlpgz:6}" == ".lp.gz" ] || [ "${line:$tmplenmpsgz:7}" == ".mps.gz" ] || [ "${line:$tmplenlpbz:7}" == ".lp.bz2" ] || [ "${line:$tmplenmpsbz:8}" == ".mps.bz2" ]; then
  echo "`date`: Running vpc code on instance ${INSTANCE_LIST}"
  echo "python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${INSTANCE_LIST} ${RESULTS_DIR} ${RUN_TYPE_STUB} \"$4\" >> ${RESULTS_DIR}/log.out 2>&1"
  python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${INSTANCE_LIST} ${RESULTS_DIR} ${RUN_TYPE_STUB} "$4" >> ${RESULTS_DIR}/log.out 2>&1
elif [ "${tmpext:0:9}" == "instances" ]; then
  echo "Using sequential mode"
  batchstub=""
  while read line; do
    # Skip empty lines
    if [ -z "$line" ]
    then
      continue
    fi

    # If line ends with "/", then it is a batch name
    len="$((${#line}-1))"
    if [ "${line:len:1}" == "/" ]; then
      batchstub="${line:0:len}"
      mkdir -p ${RESULTS_DIR}/$batchstub
      continue
    fi

    # If the line has a comma, the second half will be extra params the user defined for that instance
    inst=""
    userparams=""
    i=0
    split_lines=`split_on_commas "$line"`
    while read item; do
      if [ "$i" -eq "0" ]; then
        inst=$item
      fi
      if [ $i == 1 ]; then
        userparams=$item
      fi
      #if [ $i > 1 ]; then
      #  userparams="$userparams $item"
      #fi
      i=$((i+1))
    done <<< "$split_lines"
    
    # Add INSTANCE_DIR and .mps if no (known) extension detected
    tmplenlp="$((${#inst}-3))"
    tmplenmps="$((${#inst}-4))"
    tmplenlpgz="$((${#inst}-6))"
    tmplenmpsgz="$((${#inst}-7))"
    tmplenlpbz="$((${#inst}-7))"
    tmplenmpsbz="$((${#inst}-8))"
    if [ "${inst:$tmplenlp:3}" != ".lp" ] && [ "${inst:$tmplenmps:4}" != ".mps" ] && [ "${inst:$tmplenlpgz:6}" != ".lp.gz" ] && [ "${inst:$tmplenmpsgz:7}" != ".mps.gz" ] && [ "${inst:$tmplenlpbz:7}" != ".lp.bz2" ] && [ "${inst:$tmplenmpsbz:8}" != ".mps.bz2" ]; then
      inst="${INSTANCE_DIR}/${inst}.mps"
    fi

    #echo "Calling python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${inst} ${RESULTS_DIR} ${RUN_TYPE_STUB} \"$userparams\" >> ${RESULTS_DIR}/log.out 2>&1 &"
    #python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${inst} ${RESULTS_DIR} ${RUN_TYPE_STUB} "$userparams" >> ${RESULTS_DIR}/log.out 2>&1
    echo "Calling ${SCRIPT_DIR}/${CURRSCRIPTNAME} $inst ${RESULTS_DIR}/$batchstub ${RUN_TYPE_STUB} \"$userparams\" >> ${RESULTS_DIR}/$batchstub/log.out 2>& 1"
    ${SCRIPT_DIR}/${CURRSCRIPTNAME} ${inst} ${RESULTS_DIR}/$batchstub ${RUN_TYPE_STUB} "$userparams" >> ${RESULTS_DIR}/$batchstub/log.out 2>&1
  done < ${INSTANCE_LIST}
elif [ "${tmpext:0:5}" == "batch" ]; then
  echo "Using batch mode."
  #FSTUB=`head -n 1 ${INSTANCE_LIST}`
  FSTUB=${INSTANCE_LIST##*/}
  FSTUB=${FSTUB%.*}
  batchstub=""
  tmpfilename=""
  while read line; do
    # Skip empty lines
    if [ -z "$line" ]; then
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
        echo "Starting batch ${tmpfilename}"
        echo "Calling nohup ${SCRIPT_DIR}/${CURRSCRIPTNAME} ${tmpfilename} ${RESULTS_DIR}/$batchstub ${RUN_TYPE_STUB} \"$4\" >> ${RESULTS_DIR}/$batchstub/log.out 2>& 1 &"
        nohup ${SCRIPT_DIR}/${CURRSCRIPTNAME} ${tmpfilename} ${RESULTS_DIR}/$batchstub ${RUN_TYPE_STUB} "$4" >> ${RESULTS_DIR}/$batchstub/log.out 2>&1 &
        #echo "nohup python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${tmpfilename} ${RESULTS_DIR} ${RUN_TYPE_STUB} \"$4\" >> ${RESULTS_DIR}/log.out 2>&1 &"
        #nohup python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${tmpfilename} ${RESULTS_DIR} ${RUN_TYPE_STUB} "$4" >> ${RESULTS_DIR}/log.out 2>&1 &
      fi

      # Now we create the new batch
      batchstub="${line:0:len}"
      #echo "Creating new batch $batchstub"
      tmpfilename="/tmp/${FSTUB}.instances${batchstub}XXX"
      tmpfilename=$(mktemp -q ${tmpfilename})
      #echo "${FSTUB}" > "${tmpfilename}"
      mkdir -p ${RESULTS_DIR}/$batchstub
    else
      # Add the current line to the current batch
      echo "${line}" >> "${tmpfilename}"
    fi
  done < ${INSTANCE_LIST}

  # Now process the last batch
  if [ ! -z "${tmpfilename}" ]
  then
    echo "Starting batch ${tmpfilename}"
    echo "Calling nohup ${SCRIPT_DIR}/${CURRSCRIPTNAME} ${tmpfilename} ${RESULTS_DIR}/$batchstub ${RUN_TYPE_STUB} \"$4\" >> ${RESULTS_DIR}/$batchstub/log.out 2>& 1 &"
    nohup ${SCRIPT_DIR}/${CURRSCRIPTNAME} ${tmpfilename} ${RESULTS_DIR}/$batchstub ${RUN_TYPE_STUB} "$4" >> ${RESULTS_DIR}/$batchstub/log.out 2>&1 &
    #echo "nohup python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${tmpfilename} ${RESULTS_DIR} ${RUN_TYPE_STUB} \"$4\" >> ${RESULTS_DIR}/log.out 2>&1 &"
    #nohup python3 -u ${SCRIPT_DIR}/${SCRIPTNAME} ${tmpfilename} ${RESULTS_DIR} ${RUN_TYPE_STUB} "$4" >> ${RESULTS_DIR}/log.out 2>&1 &
  fi
else
  echo "Could not identify type of instance file given by $line"
fi
