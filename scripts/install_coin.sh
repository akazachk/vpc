#!/usr/bin/env bash
# Takes as an (optional) argument the directory where you wish to install the COIN-OR software

## User needs to define
#PROJ_DIR="~/repos/vpc"
CBC_VERSION="2.10"
#CBC_VERSION="2.9"
##CBC_REVISION="2352"
#CBC_REVISION="2376"
#CBC_URL="https://projects.coin-or.org/svn/Cbc/stable/${CBC_VERSION}"
CBC_URL="https://projects.coin-or.org/svn/Cbc/trunk"
#COIN_DIR_NAME="Cbc-${CBC_VERSION}"
COIN_DIR_NAME="Cbc-trunk"
#COIN_DIR_NAME="Cbc-${CBC_VERSION}r${CBC_REVISION}"

USE_COINBREW=1

## Check settings
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
mkdir -p $COIN_DIR

if [ ${USE_COINBREW} == "1" ]
then
  cd $COIN_DIR
  wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
  chmod u+x coinbrew
  #./coinbrew fetch Cbc:stable/2.9
  #./coinbrew fetch Cbc:stable/2.10
  ./coinbrew fetch Cbc
  #cp ${PROJ_DIR}/lib/CbcModel.* ${COIN_DIR}/Cbc/src/
  # -b: specify build directory
  # -p: install in same directory as build
  #./coinbrew build install Cbc --no-prompt -b build -p --test --with-cplex=false CXXDEFS="-DSAVE_NODE_INFO" CDEFS="-DSAVE_NODE_INFO" ADD_CXXFLAGS="-DSAVE_NODE_INFO"
  ./coinbrew build install Cbc -b build -p --no-prompt --test --with-cplex=false ADD_CXXFLAGS="-DSAVE_NODE_INFO"
  ./coinbrew build install Cbc -b buildg -p --no-prompt --test --with-cplex=false --enable-debug ADD_CXXFLAGS="-DSAVE_NODE_INFO"
else
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
  #EXTRA_ARGS=--cdefs="-DSAVE_NODE_INFO"

  if [ -z ${CBC_REVISION} ]
  then
    svn co $CBC_URL $COIN_DIR
  else
    svn co -r ${CBC_REVISION} $CBC_URL $COIN_DIR
  fi
  cd $COIN_DIR
  cp ${PROJ_DIR}/lib/CbcModel.* ${COIN_DIR}/Cbc/src/

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
fi

echo ""
echo "Done with installing COIN-OR files into $COIN_DIR"
