#!/bin/bash
# Takes as an (optional) argument the directory where you wish to install the COIN-OR software

## User needs to define
#PROJ_DIR="~/repos/vpc"
#CBC_VERSION="2.10"
CBC_VERSION="2.9"
CBC_REVISION="2352"
CBC_URL="https://projects.coin-or.org/svn/Cbc/stable/${CBC_VERSION}"
#CGL_URL="https://projects.coin-or.org/svn/Cgl/stable/0.59"
COIN_DIR_NAME="Cbc-${CBC_VERSION}"
if [ -z "$PROJ_DIR" ]
  then echo "Need to define PROJ_DIR. Exiting."
  exit
fi
if [ -z "$1" ]
then
	COIN_DIR="${PROJ_DIR}/lib/${COIN_DIR_NAME}"
else
	COIN_DIR="${1}/${COIN_DIR_NAME}"
fi

## Ignore below unless you wish to use OsiCpxSolverInterface
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
#EXTRA_ARGS=--with-cplex-incdir=$CPLEX_INC --with-cplex-lib="$CPLEX_LIB"

mkdir -p $COIN_DIR
if [ -z ${CBC_REVISION} ]
then
  svn co $CBC_URL $COIN_DIR
else
  svn co -r ${CBC_REVISION} $CBC_URL $COIN_DIR
fi
cd $COIN_DIR

mkdir -p build
cd build
../configure -C ${EXTRA_ARGS} >& last_config.txt
make
make install
cd ..

mkdir -p buildg
cd buildg
../configure -C ${EXTRA_ARGS} --enable-debug=yes >& last_config.txt
make
make install

echo ""
echo "Done with installing COIN-OR files into $COIN_DIR"
