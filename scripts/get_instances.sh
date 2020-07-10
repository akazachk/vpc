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
DATA_DIR="${VPC_DIR}/data/instances/original"

echo "VPC_DIR is set to $VPC_DIR"
echo "Instances will be downloaded to $DATA_DIR"

START_DIR=$PWD
cd $DATA_DIR

echo "Downloading MIPLIB2017"
wget http://miplib2017.zib.de/downloads/collection.zip
mkdir -p ${DATA_DIR}/miplib2017
unzip collection.zip -d ${DATA_DIR}/miplib2017
mv ${DATA_DIR}/miplib2017/mas74.mps.gz ${DATA_DIR}/miplib2017/mas074.mps.gz
mv ${DATA_DIR}/miplib2017/mas76.mps.gz ${DATA_DIR}/miplib2017/mas076.mps.gz

echo "Downloading MIPLIB2010"
wget http://miplib2010.zib.de/download/miplib2010-complete.zip
unzip miplib2010-complete.zip -d ${DATA_DIR}

echo "Downloading MIPLIB2003"
wget http://miplib2010.zib.de/miplib2003/download/miplib2003.tar
mkdir -p ${DATA_DIR}/miplib2003
tar -xvf miplib2003.tar --directory ${DATA_DIR}/miplib2003

echo "Downloading MIPLIB3"
wget http://miplib2010.zib.de/miplib3/miplib3.tar.gz
tar -xvf miplib3.tar.gz --directory ${DATA_DIR}
rm ${DATA_DIR}/miplib3/miplib.cat
rm ${DATA_DIR}/miplib3/mps_format
rm ${DATA_DIR}/miplib3/references
cd ${DATA_DIR}/miplib3
for f in *; do mv "$f" "$f.mps"; done
for f in *; do gzip "$f"; done
cd ${DATA_DIR}
cp ${DATA_DIR}/miplib2/stein27_nocard.mps.gz ${DATA_DIR}/miplib3/
cp ${DATA_DIR}/miplib2/stein45_nocard.mps.gz ${DATA_DIR}/miplib3/

#echo "Downloading MIPLIB2"
#wget http://miplib2010.zib.de/miplib2/miplib.tar.gz
#mkdir ${DATA_DIR}/miplib2
#tar -xvf miplib.tar.gz --directory ${DATA_DIR}/miplib2
#rm ${DATA_DIR}/miplib2/miplib.cat
#rm ${DATA_DIR}/miplib2/miplib/mps_format
#rm ${DATA_DIR}/miplib2/miplib/references
#cd ${DATA_DIR}/miplib2/miplib
#for f in *; do mv "$f" "$f.mps"; done
#for f in *; do gzip "$f"; done
#mv *.mps.gz ..
#cd ..
#rm -r miplib
#cd ${DATA_DIR}

echo "Downloading COR@L instances"
wget https://coral.ise.lehigh.edu/wp-content/uploads/mip-instances/instances/ALL_INSTANCE.tar
mv ALL_INSTANCE.tar coral.tar
mkdir -p ${DATA_DIR}/coral
tar -xvf coral.tar --directory ${DATA_DIR}/coral
mv coral/mcf2.mps.bz2 coral/danoint.mps.bz2

echo "Finished downloading instances to $DATA_DIR"
cd ${START_DIR}
