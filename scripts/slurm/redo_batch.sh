#!/bin/bash
#SBATCH --array=11-20,21-31
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=4G

#SBATCH --time=24:00:00
#SBATCH --array=1-16
#SBATCH --array=10-15
#SBATCH --array=21-27
#SBATCH --array=15
#SBATCH --mem-per-cpu=10G

#SBATCH --account=def-alodi
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aleksandr.kazachkov@polymtl.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

split_on_commas() {
  local IFS=,
  local WORD_LIST=($1)
  for word in "${WORD_LIST[@]}"; do
    echo "$word"
  done
}

TYPE="presolved"
MODE="bb"
CASE_NUM=`printf %03d $SLURM_ARRAY_TASK_ID`
INSTANCE_FILE=${TYPE}_redo.instances
export VPC_DIR="${REPOS_DIR}/vpc"

# Set mode if given
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

echo "Starting ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if (("$SLURM_ARRAY_TASK_ID" <= -1))
then
  echo "Running task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} in batch mode at `date`"
  FILE="${VPC_DIR}/scripts/slurm/${TYPE}${CASE_NUM}.instances"
#elif (($SLURM_ARRAY_TASK_ID >= 3)) && (($SLURM_ARRAY_TASK_ID <= 4))
#then
#  echo "Running $SLURM_ARRAY_TASK_ID in batch mode"
else
  echo "Running task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} in sequential mode at `date`"
  line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${VPC_DIR}/scripts/slurm/${INSTANCE_FILE})

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
  FILE="${VPC_DIR}/data/instances/${inst}.mps"
fi

${VPC_DIR}/scripts/run_experiments.sh $FILE ${VPC_DIR}/results/${MODE}_redo/$CASE_NUM ${MODE} $userparams

#echo "Statistics from seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
#seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
