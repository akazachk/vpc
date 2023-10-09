#!/bin/bash
#SBATCH --account=akazachkov
#SBATCH --job-name=slurm_preprocess
#SBATCH --time=03:00:00             # time limit hrs:min:sec
#SBATCH --mem-per-cpu=1G            # job memory
#SBATCH --ntasks=1                  # run a single task
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=akazachkov@ufl.edu
#SBATCH --mail-type=BEGIN,FAIL,END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=preprocess_%A-%a.log    # Standard output and error log
#SBATCH --array=1-10,329-503
#SBATCH --array=1
pwd; hostname; date

TYPE="original"
MODE="preprocess"
CASE_NUM=`printf %03d $SLURM_ARRAY_TASK_ID`
INSTANCE_FILE=${TYPE}.instances
INSTANCE_FILE=small_original.instances
export VPC_DIR="${REPOS_DIR}/vpc"
RESULTS_DIR=${VPC_DIR}/results
RESULTS_DIR="/blue/akazachkov/vpc/2022-11-01/"

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
  FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${VPC_DIR}/scripts/slurm/${INSTANCE_FILE}).mps
  FILE="${VPC_DIR}/data/instances/${FILE}"
fi

${VPC_DIR}/scripts/run_experiments.sh $FILE ${RESULTS_DIR}/${MODE}/$CASE_NUM ${MODE} $CASE_NUM

#echo "Statistics from seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
#seff ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

date
