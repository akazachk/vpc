#!/usr/bin/env bash

if [ -z $1 ]; then 
  echo "Need to give job number."
  exit 1
else
  sacct -j $1 --format=JobID,JobName,MaxRSS,Elapsed,State
fi
