#!/usr/bin/env bash

# Create batch files
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
