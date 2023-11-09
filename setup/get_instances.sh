#!/usr/bin/env bash
# Usage: get_instances.sh [data_dir] [skip_download=0/1] [all / miplib2017 / miplib2010 / miplib3 / miplib2 / coral / neos] [echo]

# Set debug mode (do not do anything; only print)
DEBUG=echo
DEBUG=

if [ -z $1 ]; then
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
  DATA_DIR="${VPC_DIR}/data/instances/original"
else
  DATA_DIR="${1}"
fi
echo "DATA_DIR is set to $DATA_DIR"

SKIP_DOWNLOAD=0
if [ ! -z $2 ]; then
  if [ $2 == 1 ]; then
    SKIP_DOWNLOAD=1
  fi
fi

SETS="all"
if [ ! -z $3 ]; then
  SETS=$3
fi

# Set debug mode (do not do anything; only print)
if [ ! -z $4 ]; then
  DEBUG=echo
fi

echo "Instances will be downloaded to $DATA_DIR"

START_DIR=$PWD
$DEBUG mkdir -p ${DATA_DIR}
#$DEBUG cd $DATA_DIR

if [ "$SETS" == "all" ] || [ "$SETS" == "miplib2017" ]; then
  INSTANCES=miplib2017
  echo ""
  echo "Processing ${INSTANCES}"
  if [ ${SKIP_DOWNLOAD} != 1 ]; then
    $DEBUG wget -O ${DATA_DIR}/${INSTANCES}.zip 'https://miplib.zib.de/downloads/collection.zip'
  fi
  $DEBUG mkdir -p ${DATA_DIR}/${INSTANCES}
  echo "Unpacking ${INSTANCES}"
  $DEBUG unzip -q ${DATA_DIR}/${INSTANCES}.zip -d ${DATA_DIR}/${INSTANCES}
  $DEBUG mv ${DATA_DIR}/${INSTANCES}/revised-submissions/*/instances/*.mps.gz ${DATA_DIR}/${INSTANCES}
  $DEBUG rm -r ${DATA_DIR}/${INSTANCES}/revised-submissions
  #mv ${DATA_DIR}/${INSTANCES}/mas74.mps.gz ${DATA_DIR}/miplib2017/mas074.mps.gz
  #mv ${DATA_DIR}/${INSTANCES}/mas76.mps.gz ${DATA_DIR}/miplib2017/mas076.mps.gz
fi

if [ "$SETS" == "all" ] || [ "$SETS" == "miplib2010" ]; then
  # In MIPLIB 2010, unzipping will create desired directory
  INSTANCES=miplib2010
  TYPE=complete
  echo ""
  echo "Processing ${INSTANCES}"
  if [ ${SKIP_DOWNLOAD} != 1 ]; then
    if [ ${TYPE} == "complete" ]; then
      $DEBUG wget -O ${DATA_DIR}/${INSTANCES}.zip 'https://miplib2010.zib.de/download/miplib2010-complete.zip'
    else
      $DEBUG wget -O ${DATA_DIR}/${INSTANCES}.zip 'https://miplib2010.zib.de/download/miplib2010-1.1.3-benchmark.zip'
    fi
  fi
  echo "Unpacking ${INSTANCES}"
  $DEBUG unzip -q ${DATA_DIR}/${INSTANCES}.zip *.mps.gz -d ${DATA_DIR}/${INSTANCES}
  if [ ${TYPE} == "complete" ]; then
    $DEBUG mv ${DATA_DIR}/${INSTANCES}/${INSTANCES}/*.mps.gz ${DATA_DIR}/${INSTANCES}
    $DEBUG rm -r ${DATA_DIR}/${INSTANCES}/miplib2010
  else
    # Did not test below
    $DEBUG mv ${DATA_DIR}/${INSTANCES}/*/instances/${INSTANCES}/*.mps.gz ${DATA_DIR}/${INSTANCES}
    $DEBUG rm -r ${DATA_DIR}/${INSTANCES}/miplib2010-*
  fi
fi

if [ "$SETS" == "all" ] || [ "$SETS" == "miplib2003" ]; then
  INSTANCES=miplib2003
  echo ""
  echo "Processing ${INSTANCES}"
  if [ ${SKIP_DOWNLOAD} != 1 ]; then
    $DEBUG wget -O ${DATA_DIR}/${INSTANCES}.tar https://miplib2010.zib.de/miplib2003/download/miplib2003.tar
  fi
  $DEBUG mkdir -p ${DATA_DIR}/${INSTANCES}
  echo "Unpacking ${INSTANCES}"
  $DEBUG tar -xf ${DATA_DIR}/${INSTANCES}.tar --directory ${DATA_DIR}/${INSTANCES}
fi

if [ "$SETS" == "all" ] || [ "$SETS" == "miplib3" ]; then
  INSTANCES=miplib3
  echo ""
  echo "Processing ${INSTANCES}"
  if [ ${SKIP_DOWNLOAD} != 1 ]; then
    $DEBUG wget -O ${DATA_DIR}/${INSTANCES}.tar.gz https://miplib2010.zib.de/miplib3/miplib3.tar.gz
  fi
  echo "Unpacking ${INSTANCES}"
  $DEBUG tar -xf ${DATA_DIR}/${INSTANCES}.tar.gz --directory ${DATA_DIR}
  $DEBUG rm ${DATA_DIR}/${INSTANCES}/miplib.cat
  $DEBUG rm ${DATA_DIR}/${INSTANCES}/mps_format
  $DEBUG rm ${DATA_DIR}/${INSTANCES}/references
  CURR_DIR=$PWD
  cd ${DATA_DIR}/${INSTANCES}
  for f in *; do $DEBUG mv "$f" "$f.mps"; done
  for f in *; do $DEBUG gzip "$f"; done
  cd ${CURR_DIR}
fi

if [ "$SETS" == "all" ] || [ "$SETS" == "miplib2" ]; then
  INSTANCES=miplib2
  echo ""
  echo "Processing ${INSTANCES}"
  if [ ${SKIP_DOWNLOAD} != 1 ]; then
    $DEBUG wget -O ${DATA_DIR}/${INSTANCES}.tar.gz 'https://miplib2010.zib.de/miplib2/miplib.tar.gz'
  fi
  $DEBUG mkdir ${DATA_DIR}/${INSTANCES}
  echo "Unpacking ${INSTANCES}"
  $DEBUG tar -xf ${DATA_DIR}/${INSTANCES}.tar.gz --directory ${DATA_DIR}/${INSTANCES}
  $DEBUG mv ${DATA_DIR}/${INSTANCES}/miplib/* ${DATA_DIR}/${INSTANCES}
  $DEBUG rm ${DATA_DIR}/${INSTANCES}/miplib.cat
  $DEBUG rm ${DATA_DIR}/${INSTANCES}/mps_format
  $DEBUG rm ${DATA_DIR}/${INSTANCES}/references
  CURR_DIR=$PWD
  cd ${DATA_DIR}/${INSTANCES}
  $DEBUG rm -r miplib
  for f in *; do $DEBUG mv "$f" "$f.mps"; done
  for f in *; do $DEBUG gzip "$f"; done
  #$DEBUG mv *.mps.gz ..
  #$DEBUG cd ..
  cd ${CURR_DIR}
  #cp ${DATA_DIR}/miplib2/stein27_nocard.mps.gz ${DATA_DIR}/miplib3/
  #cp ${DATA_DIR}/miplib2/stein45_nocard.mps.gz ${DATA_DIR}/miplib3/
fi

if [ "$SETS" == "all" ] || [ "$SETS" == "coral" ]; then
  INSTANCES=coral
  echo ""
  echo "Processing ${INSTANCES}"
  if [ ${SKIP_DOWNLOAD} != 1 ]; then
    $DEBUG wget -O ${DATA_DIR}/${INSTANCES}.tar https://coral.ise.lehigh.edu/wp-content/uploads/mip-instances/instances/ALL_INSTANCE.tar
  fi
  $DEBUG mkdir -p ${DATA_DIR}/${INSTANCES}
  echo "Unpacking ${INSTANCES}"
  $DEBUG tar -xf ${DATA_DIR}/${INSTANCES}.tar --directory ${DATA_DIR}/${INSTANCES}
  $DEBUG mv ${DATA_DIR}/${INSTANCES}/mcf2.mps.bz2 ${DATA_DIR}/${INSTANCES}/danoint.mps.bz2
fi

if [ "$SETS" == "all" ] || [ "$SETS" == "neos" ]; then
  INSTANCES=neos
  echo ""
  echo "Processing ${INSTANCES}"
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
  $DEBUG cp -r ${VPC_DIR}/data/instances/original/neos ${DATA_DIR}/${INSTANCES}
fi

echo "Finished downloading instances to $DATA_DIR"
cd ${START_DIR}
