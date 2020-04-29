#!/usr/bin/env bash

if [ -z "$PROJ_DIR" ]
  then echo "Need to define PROJ_DIR. Exiting."
  exit
fi
DATA_DIR="${PROJ_DIR}/data/instances"

wget http://miplib2017.zib.de/downloads/collection.zip
mkdir -p ${DATA_DIR}/original/miplib2017
unzip collection.zip -d ${DATA_DIR}/original/miplib2017

wget http://miplib2010.zib.de/download/miplib2010-complete.zip
unzip miplib2010-complete.zip -d ${DATA_DIR}/original

wget http://miplib2010.zib.de/miplib2003/download/miplib2003.tar
mkdir -p ${DATA_DIR}/original/miplib2003
tar -xvf miplib2003.tar --directory ${DATA_DIR}/original/miplib2003

wget http://miplib2010.zib.de/miplib3/miplib3.tar.gz
tar -xvf miplib3.tar.gz --directory ${DATA_DIR}/original
rm ${DATA_DIR}/original/miplib3/miplib.cat
rm ${DATA_DIR}/original/miplib3/mps_format
rm ${DATA_DIR}/original/miplib3/references
cd ${DATA_DIR}/original/miplib3
for f in *; do mv "$f" "$f.mps"; done
for f in *; do gzip "$f"; done
cd ${PROJ_DIR}

wget http://miplib2010.zib.de/miplib2/miplib.tar.gz
mkdir ${DATA_DIR}/original/miplib2
tar -xvf miplib.tar.gz --directory ${DATA_DIR}/original/miplib2
rm ${DATA_DIR}/original/miplib2/miplib.cat
rm ${DATA_DIR}/original/miplib2/miplib/mps_format
rm ${DATA_DIR}/original/miplib2/miplib/references
cd ${DATA_DIR}/original/miplib2/miplib
for f in *; do mv "$f" "$f.mps"; done
for f in *; do gzip "$f"; done
mv *.mps.gz ..
cd ..
rm -r miplib
cd ${PROJ_DIR}
