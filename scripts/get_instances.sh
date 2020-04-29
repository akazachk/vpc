#!/usr/bin/env bash

wget http://miplib2017.zib.de/downloads/collection.zip
mkdir -p original/miplib2017
unzip collection.zip -d original/miplib2017

wget http://miplib2010.zib.de/download/miplib2010-complete.zip
unzip miplib2010-complete.zip -d original

wget http://miplib2010.zib.de/miplib2003/download/miplib2003.tar
mkdir -p original/miplib2003
tar -xvf miplib2003.tar --directory original/miplib2003

wget http://miplib2010.zib.de/miplib3/miplib3.tar.gz
tar -xvf miplib3.tar.gz --directory original
rm original/miplib3/miplib.cat
rm original/miplib3/mps_format
rm original/miplib3/references
cd original/miplib3
for f in *; do mv "$f" "$f.mps"; done
for f in *; do gzip "$f"; done
cd ../..

wget http://miplib2010.zib.de/miplib2/miplib.tar.gz
mkdir original/miplib2
tar -xvf miplib.tar.gz --directory original/miplib2
rm original/miplib2/miplib.cat
rm original/miplib2/miplib/mps_format
rm original/miplib2/miplib/references
cd original/miplib2/miplib
for f in *; do mv "$f" "$f.mps"; done
for f in *; do gzip "$f"; done
mv *.mps.gz ..
cd ..
rm -r miplib
cd ../..
