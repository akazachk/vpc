#!/usr/bin/env bash

# Create batch files
TYPE="presolved"
if [ $TYPE = "original" ]; then
  echo "original!"
  SRC=original.instances
  EXT="instances"
  head -n 160 $SRC > 1.${EXT}
  tail -n +161 $SRC | head -n 25 > 2.${EXT}
  tail -n +186 $SRC | head -n 15 > 3.${EXT}
  tail -n +201 $SRC | head -n 6 > 4.${EXT}
  tail -n +207 $SRC | head -n 5 > 5.${EXT}
  tail -n +212 $SRC | head -n 5 > 6.${EXT}
  tail -n +217 $SRC | head -n 3 > 7.${EXT}
  tail -n +220 $SRC | head -n 3 > 8.${EXT}
  tail -n +223 $SRC | head -n 3 > 9.${EXT}
  exit 1
elif [ $TYPE = "presolved" ]; then
  echo "presolved!"
  SRC=${TYPE}.instances
  EXT="instances"
  head -n 125 $SRC > 1.${EXT}
  tail -n +126 $SRC | head -n 15 > 2.${EXT}
  tail -n +141 $SRC | head -n 10 > 3.${EXT}
  tail -n +151 $SRC | head -n 6 > 4.${EXT}
  tail -n +157 $SRC | head -n 6 > 5.${EXT}
  tail -n +164 $SRC | head -n 5 > 6.${EXT}
  tail -n +168 $SRC | head -n 4 > 7.${EXT}
  tail -n +172 $SRC | head -n 3 > 8.${EXT}
  tail -n +175 $SRC | head -n 3 > 9.${EXT}
  tail -n +178 $SRC | head -n 3 > 10.${EXT}
  exit 1
else
  echo "Unrecognized type: $TYPE"
  exit 1
fi

