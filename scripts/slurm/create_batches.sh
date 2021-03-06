#!/usr/bin/env bash

# Create batch files
TYPE="presolved"
#TYPE="presolved"

if [ ! -z $1 ]
then TYPE=$1
fi

if [ $TYPE = "original" ]; then
  echo "original!"
  SRC=original.instances
  EXT="instances"
  head -n 167 $SRC > ${TYPE}001.${EXT}
  tail -n +168 $SRC | head -n 78 > ${TYPE}011.${EXT}
  tail -n +246 $SRC | head -n 14 > ${TYPE}002.${EXT}
  tail -n +260 $SRC | head -n 16 > ${TYPE}012.${EXT}
  tail -n +276 $SRC | head -n 4 > ${TYPE}003.${EXT}
  tail -n +280 $SRC | head -n 6 > ${TYPE}013.${EXT}
  tail -n +286 $SRC | head -n 5 > ${TYPE}014.${EXT}
  tail -n +291 $SRC | head -n 6 > ${TYPE}004.${EXT}
  tail -n +297 $SRC | head -n 4 > ${TYPE}015.${EXT}
  tail -n +301 $SRC | head -n 7 > ${TYPE}005.${EXT}
  tail -n +308 $SRC | head -n 6 > ${TYPE}006.${EXT}
  tail -n +314 $SRC | head -n 3 > ${TYPE}007.${EXT}
  tail -n +317 $SRC | head -n 2 > ${TYPE}016.${EXT}
  tail -n +319 $SRC | head -n 4 > ${TYPE}008.${EXT}
  tail -n +323 $SRC | head -n 3 > ${TYPE}009.${EXT}
  tail -n +326 $SRC | head -n 3 > ${TYPE}010.${EXT}
  exit 1
elif [ $TYPE = "presolved" ]; then
  echo "presolved!"
  SRC=${TYPE}.instances
  EXT="instances"
  head -n 89 $SRC > ${TYPE}001.${EXT}
  tail -n +90 $SRC | head -n 35 > ${TYPE}011.${EXT}
  tail -n +126 $SRC | head -n 15 > ${TYPE}002.${EXT}
  tail -n +141 $SRC | head -n 10 > ${TYPE}003.${EXT}
  tail -n +151 $SRC | head -n 6 > ${TYPE}004.${EXT}
  tail -n +157 $SRC | head -n 6 > ${TYPE}005.${EXT}
  tail -n +164 $SRC | head -n 5 > ${TYPE}006.${EXT}
  tail -n +168 $SRC | head -n 4 > ${TYPE}007.${EXT}
  tail -n +172 $SRC | head -n 3 > ${TYPE}008.${EXT}
  tail -n +175 $SRC | head -n 3 > ${TYPE}009.${EXT}
  tail -n +178 $SRC | head -n 2 > ${TYPE}010.${EXT}
  exit 1
else
  echo "Unrecognized type: $TYPE"
  exit 1
fi

