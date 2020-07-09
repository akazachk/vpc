#!/bin/bash
#SBATCH --array=1-12,229-489
#SBATCH --time=30:00:00
#SBATCH --account=def-alodi
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aleksandr.kazachkov@polymtl.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

MODE="bb"
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
INSTANCE_FILE=original.instances

if (("$SLURM_ARRAY_TASK_ID" <= 12))
then
  echo "Running task $SLURM_ARRAY_TASK_ID in batch mode at `date`"
  FILE=${REPOS_DIR}/vpc/scripts/slurm/${SLURM_ARRAY_TASK_ID}.batch
#elif (($SLURM_ARRAY_TASK_ID >= 3)) && (($SLURM_ARRAY_TASK_ID <= 4))
#then
#  echo "Running $SLURM_ARRAY_TASK_ID in batch mode"
else
  echo "Running task $SLURM_ARRAY_TASK_ID in sequential mode at `date`"
  FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${REPOS_DIR}/vpc/scripts/slurm/${INSTANCE_FILE})
fi

${REPOS_DIR}/vpc/scripts/run_experiments.sh $FILE ${REPOS_DIR}/vpc/results/${MODE}/$CASE_NUM ${MODE} $CASE_NUM

echo "Statistics from seff $SLURM_ARRAY_JOB_ID_$SLURM_ARRAY_TASK_ID"
seff $SLURM_ARRAY_JOB_ID_$SLURM_ARRAY_TASK_ID
