#!/usr/bin/env bash

# Create batch files
SRC=original.instances
head -n 160 $SRC > 1.batch
tail -n +161 $SRC | head -n 25 > 2.batch
tail -n +186 $SRC | head -n 15 > 3.batch
tail -n +201 $SRC | head -n 5 > 4.batch
tail -n +206 $SRC | head -n 5 > 5.batch
tail -n +211 $SRC | head -n 3 > 6.batch
tail -n +214 $SRC | head -n 3 > 7.batch
tail -n +217 $SRC | head -n 3 > 8.batch
tail -n +220 $SRC | head -n 3 > 9.batch
tail -n +223 $SRC | head -n 2 > 10.batch
tail -n +225 $SRC | head -n 2 > 11.batch
tail -n +227 $SRC | head -n 2 > 12.batch
