#!/usr/bin/env bash
# Takes as an (optional) argument the directory where you wish to install the COIN-OR software

## User needs to define
#VPC_DIR="~/repos/vpc"
#CBC_VERSION="2.10"
#CBC_VERSION="2.9"
CBC_VERSION="trunk"
##CBC_REVISION="2352"
#CBC_REVISION="2376"
if [ ${CBC_VERSION} == "trunk" ]
then
  CBC_URL="https://projects.coin-or.org/svn/Cbc/trunk"
else
  CBC_URL="https://projects.coin-or.org/svn/Cbc/stable/${CBC_VERSION}"
fi
COIN_DIR_NAME="Cbc-${CBC_VERSION}"
#COIN_DIR_NAME="Cbc-${CBC_VERSION}r${CBC_REVISION}"

USE_COINBREW=1

## Check settings
if [ -z "$CPLEX_HOME" ]
then
  echo "[optional] Please locate CPLEX if you would like to use it (e.g., /Applications/CPLEX_Studio1210/cplex/ on mac), or press enter if you do not wish to use CPLEX:"
  read CPLEX_HOME
fi
if [ ! -z "$CPLEX_HOME" ]
then
  echo "CPLEX_HOME set to $CPLEX_HOME"
fi

if [ -z "$GUROBI_HOME" ]
then
  echo "[optional] Please locate Gurobi if you would like to use it (e.g., /Library/gurobi_902/mac64 on mac), or press enter if you do not wish to use Gurobi:"
  read GUROBI_HOME
fi
if [ ! -z "$GUROBI_HOME" ]
then
  echo "GUROBI_HOME set to $GUROBI_HOME"
fi

if [ -z "$1" ]
then
  if [ -z "$VPC_DIR" ]
  then
    if [ ! -z "${REPOS_DIR}" ]
    then
      echo "COIN-OR files will be installed into VPC_DIR/lib/${COIN_DIR_NAME}. Please define VPC_DIR (the root vpc dir, possibly ${REPOS_DIR}/vpc):"
    else
      echo "COIN-OR files will be installed into VPC_DIR/lib/${COIN_DIR_NAME}. Please define VPC_DIR (the root vpc dir):"
    fi
    read VPC_DIR
    echo "Set VPC_DIR=$VPC_DIR"
    if [ -z "$VPC_DIR" ]
      then echo "Need to define VPC_DIR. Exiting."
      exit
    fi
  fi
  COIN_DIR="${VPC_DIR}/lib/${COIN_DIR_NAME}"
else
  COIN_DIR="${1}/${COIN_DIR_NAME}"
fi
echo "COIN-OR files will be installed into ${COIN_DIR}"
echo ""
#exit
mkdir -p $COIN_DIR

if [ ${USE_COINBREW} == "1" ]
then
  cd $COIN_DIR
  COINBREW_LINK="https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew"
  if which wget > /dev/null ; then
    wget -N $COINBREW_LINK
  elif which curl > /dev/null ; then
    curl $COINBREW_LINK -o coinbrew
  else
    error "Cannot download, neither wget nor curl is available."
  fi
  chmod u+x coinbrew
  if [ ${CBC_VERSION} == "trunk" ]
  then
    ./coinbrew fetch Cbc@master
  else
    ./coinbrew fetch Cbc@${CBC_VERSION}
  fi
  # -b: specify build directory
  # -p: install in same directory as build
  #./coinbrew build install Cbc -b build -p --no-prompt --test --with-cplex=false ADD_CXXFLAGS="-DSAVE_NODE_INFO"
  #./coinbrew build install Cbc -b buildg -p --no-prompt --test --with-cplex=false --enable-debug ADD_CXXFLAGS="-DSAVE_NODE_INFO"
  ./coinbrew build Cbc -b build -p build --no-prompt ADD_CXXFLAGS="-DSAVE_NODE_INFO"
  ./coinbrew build Cbc -b buildg -p buildg --no-prompt --enable-debug ADD_CXXFLAGS="-DSAVE_NODE_INFO"
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
