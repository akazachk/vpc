#!/bin/bash
# Takes as an (optional) argument the directory where you wish to install the COIN-OR software

#CGL_URL="https://projects.coin-or.org/svn/Cgl/stable/0.59"
CBC_URL="https://projects.coin-or.org/svn/Cbc/stable/2.9"
COIN_DIR_NAME="Cbc-2.9"
if [ -z "$1" ]
then
	COIN_DIR="lib/${COIN_DIR_NAME}"
else
	COIN_DIR="${1}/lib/${COIN_DIR_NAME}"
fi
#UNAME=`uname`
#if [ "$UNAME" = "Darwin" ]
#then
#	CPLEX_ARCH="x86-64_osx"
#	CPLEX_DIR="/Applications/CPLEX_Studio128"
#else
#	CPLEX_ARCH=x86-64_linux
#	CPLEX_DIR=/home/ibm/cplex-studio/12.8.0.0
#fi
#CPLEX_INC="$CPLEX_DIR/cplex/include/ilcplex"
#CPLEX_LIB_DIR="$CPLEX_DIR/cplex/lib/$CPLEX_ARCH/static_pic"
#CPLEX_LIB="-L$CPLEX_LIB_DIR -lcplex -lm -ldl -lpthread"

mkdir -p $COIN_DIR
svn co $CBC_URL $COIN_DIR
cd $COIN_DIR
mkdir -p build
cd build
#../configure -C --with-cplex-incdir=$CPLEX_INC --with-cplex-lib="$CPLEX_LIB" >& last_config.txt
../configure -C >& last_config.txt
make
make install

cd ..
mkdir -p buildg
cd buildg
#../configure -C --with-cplex-incdir=$CPLEX_INC --with-cplex-lib="$CPLEX_LIB" --enable-debug=yes >& last_config.txt
../configure -C --enable-debug=yes >& last_config.txt
make
make install

echo ""
echo "Done with installing COIN-OR files into $COIN_DIR"
