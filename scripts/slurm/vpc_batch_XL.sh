#!/bin/bash
#SBATCH --array=354-488
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=1G

#SBATCH --array=299-300
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G

#SBATCH --account=def-alodi
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aleksandr.kazachkov@polymtl.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

TYPE="presolved"
MODE="bb"
CASE_NUM=`printf %03d $SLURM_ARRAY_TASK_ID`
INSTANCE_FILE=${TYPE}.instances
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
  FILE="${VPC_DIR}/scripts/slurm/${SLURM_ARRAY_TASK_ID}.instances"
#elif (($SLURM_ARRAY_TASK_ID >= 3)) && (($SLURM_ARRAY_TASK_ID <= 4))
#then
#  echo "Running $SLURM_ARRAY_TASK_ID in batch mode"
else
  echo "Running task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} in sequential mode at `date`"
  FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${VPC_DIR}/scripts/slurm/${INSTANCE_FILE}).mps
  FILE="${VPC_DIR}/data/instances/${FILE}"
fi

${VPC_DIR}/scripts/run_experiments.sh $FILE ${VPC_DIR}/results/${MODE}/$CASE_NUM ${MODE}

#echo "Statistics from seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
#seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
