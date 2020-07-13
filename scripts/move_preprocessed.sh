#!/usr/bin/env bash

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

SOLVER_TYPE="gurobi"
SRC_DIR="${VPC_DIR}/data/instances/original"
DEST_DIR="${VPC_DIR}/data/instances/presolved_${SOLVER_TYPE}"

STUB="miplib2"
echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
mkdir -p ${DEST_DIR}/${STUB}
mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}

STUB="miplib3"
echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
mkdir -p ${DEST_DIR}/${STUB}
mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}

STUB="miplib2003"
echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
mkdir -p ${DEST_DIR}/${STUB}
mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}

STUB="miplib2010"
echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
mkdir -p ${DEST_DIR}/${STUB}
mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}

STUB="miplib2017"
echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
mkdir -p ${DEST_DIR}/${STUB}
mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}

STUB="coral"
echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
mkdir -p ${DEST_DIR}/${STUB}
mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}

STUB="neos"
echo "Copying preprocessed instances from ${SRC_DIR}/${STUB} to ${DEST_DIR}/${STUB}"
mkdir -p ${DEST_DIR}/${STUB}
mv ${SRC_DIR}/${STUB}/*_presolved.mps* ${DEST_DIR}/${STUB}
