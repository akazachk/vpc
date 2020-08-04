#!/bin/bash
#SBATCH --array=300-303
#SBATCH --time=24:00:00

#SBATCH --array=334,337,340,343,346,348,351-354,357-359,365-367,370-372,378-381,383,389-395,397,398,400-409,413,415,418,420-423,425,427-433,436,439,446,447,451,453,454,456,458,459,461-463,466,470-475,477,479-483,485,486,488,492-494,497,500,502,503
#SBATCH --time=03:00:00

#SBATCH --mem-per-cpu=10G

#SBATCH --account=def-alodi
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aleksandr.kazachkov@polymtl.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

TYPE="original"
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
