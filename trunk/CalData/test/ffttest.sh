#! /bin/bash

experiment=''
if [ $# -gt 0 ]; then
  experiment=$1
fi

if [ x$experiment = x ]; then
  echo "Usage: ffttest.sh <experiment>"
  exit 1
fi

if [ $experiment = microboone ];then
  experiment=uboone
fi

fcl1=simwire_${experiment}.fcl
if [ ! -f $fcl1 ]; then
  echo "$fcl1 not found."
  exit 1
fi

fcl2=ffttest_${experiment}.fcl
if [ ! -f $fcl2 ]; then
  echo "$fcl2 not found."
  exit 1
fi

lar -c $fcl1
lar -c $fcl2
