#!/usr/bin/env bash
# Takes one argument: gurobi or cplex (to specify type of preprocessed instances we are saving)

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
    exit
  fi
fi

if [ -z "$1" ]; then
  SOLVER_TYPE="gurobi"
else
  SOLVER_TYPE="$1"
fi

SRC_DIR="${VPC_DIR}/data/instances/original"
DEST_DIR="${VPC_DIR}/data/instances/presolved_${SOLVER_TYPE}"
SOL_DIR="${VPC_DIR}/data/solutions"
STUBS=("miplib2" "miplib3" "miplib2003" "miplib2010" "miplib2017" "coral" "neos")

for STUB in ${STUBS[*]}; do
  echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
  mkdir -p ${DEST_DIR}/${STUB}
  mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}
  mkdir -p ${SOL_DIR}/${STUB}
  if [ $SOLVER_TYPE = "cplex" ]; then
    mv ${SRC_DIR}/${STUB}/*_presolved.pre ${SOL_DIR}/${STUB}
  fi
  if [ $SOLVER_TYPE = "gurobi" ]; then
    mv ${SRC_DIR}/${STUB}/*_presolved_gurobi.mst* ${SOL_DIR}/${STUB}
  fi
done
