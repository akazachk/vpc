#/bin/bash

if [ -z ${PROJ_DIR} ]
  then echo "Need to define PROJ_DIR. Exiting."
  exit
fi

${PROJ_DIR}/Release/vpc -f ${PROJ_DIR}/test/bm23.mps -d 2 --optfile=${PROJ_DIR}/data/ip_obj.csv $1 $2 $3 $4 $5 $6 $7 $8
