#!/usr/bin/env bash
# Example usage:
#   move_preprocessed.sh ${VPC_DIR}/data/instances/original gurobi
#
# Takes three optional arguments:
#   First argument is SRC_DIR location of original files (default is $VPC_DIR/data/instances/original),
#     which also determines the destination of the presolved files (SRC_DIR/../presolved_${SOLVER_TYPE})
#   Second argument is gurobi or cplex to specify type (SOLVER_TYPE) of preprocessed instances we are saving (default is gurobi)
#   Third argument is SOL_DIR, where to copy solutions to (default is $VPC_DIR/data/solutions)


if [ -z "$1" ]; then
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
  SRC_DIR="${VPC_DIR}/data/instances/original"
else
  SRC_DIR="$1"
fi

if [ -z "$2" ]; then
  SOLVER_TYPE="gurobi"
else
  SOLVER_TYPE="$2"
fi

if [ -z "$3" ]; then
  SOL_DIR="${SRC_DIR}/../../solutions"
else
  SOL_DIR="$3"
fi

DEST_DIR="${SRC_DIR}/../presolved_${SOLVER_TYPE}"
STUBS=("miplib2" "miplib3" "miplib2003" "miplib2010" "miplib2017" "coral" "neos")

for STUB in ${STUBS[*]}; do
  echo "Moving preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
  mkdir -p ${DEST_DIR}/${STUB}
  mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}
  mkdir -p ${SOL_DIR}/${STUB}
  if [ $SOLVER_TYPE = "cplex" ]; then
    EXT="_presolved.pre"
    echo "Moving ${SOLVER_TYPE} solutions (*${EXT}) from ${SRC_DIR}/${STUB} to ${SOL_DIR}/${STUB}"
    mv ${SRC_DIR}/${STUB}/*${EXT} ${SOL_DIR}/${STUB}
  fi
  if [ $SOLVER_TYPE = "gurobi" ]; then
    EXT="_presolved_gurobi.mst"
    echo "Moving ${SOLVER_TYPE} solutions (*${EXT}*) from ${SRC_DIR}/${STUB} to ${SOL_DIR}/${STUB}"
    mv ${SRC_DIR}/${STUB}/*${EXT}* ${SOL_DIR}/${STUB}
  fi
done
