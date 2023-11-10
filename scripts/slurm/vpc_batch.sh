#!/bin/bash
#SBATCH --job-name=vpc-preprocess
#SBATCH --output=slurm_%x_%A_%a.log # standard output and error log

#SBATCH --ntasks=1                  # run a single task
#SBATCH --cpus-per-task=1

#SBATCH --mail-type=BEGIN           # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=akazachkov@ufl.edu

#SBATCH --account=akazachkov
#SBATCH --qos=akazachkov

#SBATCH --time=05:00:00             # time limit hrs:min:sec
#SBATCH --mem-per-cpu=4G            # job memory
##SBATCH --array=1-659
##SBATCH --array=1-30,32-44,46-73,75-163,165-268,270-294,296-298,301-307,309-365,367-595,597-599,601-615,617-659

#SBATCH --mem-per-cpu=8G            # job memory
#SBATCH --array=45,164,295,616

#SBATCH --mem-per-cpu=32G            # job memory
#SBATCH --array=31,74,299,300,600

#SBATCH --mem-per-cpu=8G            # job memory
#SBATCH --array=269,308,366,596




#########################
## To run this script, call (for example)
##     sbatch vpc_batch.sh bb presolved
## See arguments below
echo "=== START SLURM SCRIPT MESSAGES ==="
pwd; hostname; date


#########################
## Arguments
# Argument 1: MODE sets preset options to use; can take values "test", "bb", "bb0", or "preprocess"
MODE="bb"
# Argument 2: TYPE determines whether the "original" or "presolved" instance batch file will be used; can take values "test", "original", or "presolved"
TYPE="presolved"

# Set mode, if given
if [ ! -z $1 ]
then
  if [ $1 == "test" ] || [ $1 == "bb" ] || [ $1 == "bb0" ] || [ $1 == "preprocess" ]
  then
    MODE=$1
  else
    echo "Unrecognized mode: $1. Exiting."
    exit
  fi
fi

# Set type of instances to use, if given
if [ ! -z $2 ]; then
  if [ $2 == "test" ] || [ $2 == "original" ] || [ $2 == "presolved" ]; then
    TYPE=$2
  else
    echo "Unrecognized type: $2. Exiting."
    exit
  fi
fi


#########################
## Constants
if [ -z $REPOS_DIR ]; then
  echo "REPOS_DIR environment variable needs to be defined as the parent directory where the project is located."
  exit
fi
export PROJ_DIR="${REPOS_DIR}/vpc"
CASE_NUM=`printf %03d $SLURM_ARRAY_TASK_ID`
# LOCAL_DIR is based on machine
LOCAL_DIR=/blue/akazachkov/$USER
# INSTANCE_DIR is the parent directory where all instances are located
INSTANCE_DIR=${PROJ_DIR}/data/instances
# HiPerGator INSTANCE_DIR
INSTANCE_DIR=${LOCAL_DIR}/instances/vpc
# Where scripts are located
SCRIPT_DIR=${PROJ_DIR}/scripts
# Where instance list is located
INSTANCE_FILE="${PROJ_DIR}/experiments/${TYPE}.test"
# Where to put log files
RESULTS_DIR=${LOCAL_DIR}/results
STUB=`date +%F`
OUT_DIR=${RESULTS_DIR}/$STUB/${MODE}/${CASE_NUM}


#########################
## Prepare run command
echo "Starting ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if (("$SLURM_ARRAY_TASK_ID" <= -1))
then
  # Do a few instances sequentially, if they are fast
  # Put these in a file named ${TYPE}001.instances, etc., and modify the -1 above
  echo "Running task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} in batch mode at `date`"
  FILE="${SCRIPT_DIR}/${TYPE}${CASE_NUM}.instances"
else
  # Run one instance at a time, in parallel with other tasks in the array
  # Read one line from the provided instance file, based on the array task id
  echo "Running task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} in sequential mode at `date`"
  FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INSTANCE_FILE}).mps
  FILE="${INSTANCE_DIR}/${FILE}"
fi

echo "Calling ${SCRIPT_DIR}/run_experiments.sh $FILE ${OUT_DIR} ${MODE}"
echo "=== END SLURM SCRIPT MESSAGES ==="
echo ""


#########################
## RUN COMMAND HERE
${SCRIPT_DIR}/run_experiments.sh $FILE ${OUT_DIR} ${MODE}


#########################
## Wrap up
echo ""
#echo "Statistics from seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
#seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
