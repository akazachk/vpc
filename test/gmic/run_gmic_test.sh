#!/usr/bin/env bash

if [ -z "$VPC_DIR" ]
then 
  if [ -z "${REPOS_DIR}" ]
  then
    echo "Please define VPC_DIR (the root vpc dir):"
  else
    echo "Please define VPC_DIR (the root vpc dir, possibly ${REPOS_DIR}/vpc):"
  fi
  read VPC_DIR
  echo "Set VPC_DIR=$VPC_DIR"
  if [ -z "$VPC_DIR" ]
    then echo "Need to define VPC_DIR. Exiting."
    exit
  fi
fi

# Cbc: --bb_runs=1 --bb_mode=11 --bb_strategy=528
# CPLEX: --bb_runs=1 --bb_mode=11 --bb_strategy=532
# Gurobi: --bb_runs=1 --bb_mode=11 --bb_strategy=536
${VPC_DIR}/Debug/vpc -f ${VPC_DIR}/test/gmic/23588_presolved.mps --optfile=${VPC_DIR}/data/ip_obj.csv --logfile=${VPC_DIR}/test/gmic/log-gmic.csv --disj_terms=0 --temp=16 --rounds=100 --timelimit=1000 --gomory=1 -v0 $1 $2 $3 $4 $5 $6 $7 $8
