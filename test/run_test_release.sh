#/bin/bash

if [ -z ${PROJ_DIR} ]
  then echo "Need to define PROJ_DIR. Exiting."
  exit
fi

${PROJ_DIR}/Release/vpc -f ${PROJ_DIR}/test/bm23.mps -d 2 --optfile=${PROJ_DIR}/data/ip_obj.csv --rounds=1 -t 3600 --bb_runs=1 --bb_mode=10 --use_all_ones=1 --use_iter_bilinear=1 --use_disj_lb=1 --use_tight_points=0 --use_tight_rays=0 --use_unit_vectors=0 --logfile=log.csv $1 $2 $3 $4 $5 $6 $7 $8
#${PROJ_DIR}/Release/vpc -f ${PROJ_DIR}/test/bm23.mps -d 2 --optfile=${PROJ_DIR}/data/ip_obj.csv $1 $2 $3 $4 $5 $6 $7 $8
