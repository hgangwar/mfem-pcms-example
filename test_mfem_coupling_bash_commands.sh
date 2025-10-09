#!/usr/bin/env bash
# This script is used to run the MFEM coupling test with mpirun

# load the sources
#source ../loads.sh
PROGRAM_NAME=$1
echo $PROGRAM_NAME

#gdb --args 
mpirun -np 3 $PROGRAM_NAME 
OUTPUT=$?

if [ $OUTPUT -eq 0 ]
then
  echo "Successful!"
else
  echo "Failed!"
fi
