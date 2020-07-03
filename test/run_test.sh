#!/usr/bin/env bash

if [ -z "$VPC_DIR" ]
then 
  if [ -z "${REPOS_DIR}" ]
  then
    echo "Please define VPC_DIR (the root vpc dir, possibly ${REPOS_DIR}/vpc):"
  else
    echo "Please define VPC_DIR (the root vpc dir):"
  fi
  read VPC_DIR
  echo "Set VPC_DIR=$VPC_DIR"
  if [ -z "$VPC_DIR" ]
    then echo "Need to define VPC_DIR. Exiting."
    exit
  fi
fi

${VPC_DIR}/Debug/vpc -f ${VPC_DIR}/test/bm23.mps -d 2 --optfile=${VPC_DIR}/data/ip_obj.csv $1 $2 $3 $4 $5 $6 $7 $8
