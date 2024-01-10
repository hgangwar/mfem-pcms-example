#!/usr/bin/env bash
# This script is used to run the MFEM coupling test with mpirun

# load the sources
source ../loads.sh
PROGRAM_NAME=$1
#echo $PROGRAM_NAME

mpirun -np 3 $PROGRAM_NAME
OUTPUT=$?

if [ $OUTPUT -eq 0 ]
then
  exit 0
else
  exit 1
fi
