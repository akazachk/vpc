#!/bin/bash
#SBATCH --array=1-2,4
#SBATCH --time=00:05:00
#SBATCH --account=def-alodi
#SBATCH --mem-per-cpu=1M
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aleksandr.kazachkov@polymtl.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

MODE="test"
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

CASE_NUM=`printf %03d $SLURM_ARRAY_TASK_ID`
INSTANCE_FILE="test.instances"
export VPC_DIR="${REPOS_DIR}/vpc"

echo "Starting ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
if (("$SLURM_ARRAY_TASK_ID" <= 2))
then
  echo "Running task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} in batch mode at `date`"
  FILE="${VPC_DIR}/scripts/slurm/test${SLURM_ARRAY_TASK_ID}.instances"
else
  echo "Running task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} in sequential mode at `date`"
  FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${VPC_DIR}/scripts/slurm/${INSTANCE_FILE}).mps
  FILE="${VPC_DIR}/data/instances/${FILE}"
fi

${VPC_DIR}/scripts/run_experiments.sh $FILE ${VPC_DIR}/results/${MODE}/$CASE_NUM ${MODE} $CASE_NUM

echo "Statistics from seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
