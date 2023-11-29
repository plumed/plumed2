#! /usr/bin/env bash

if [[ -z $PLUMED_KERNEL ]]; then
  echo "$(basename $0) can work only if \"PLUMED_KERNEL\" is defined"
  echo "either via module load or sourceme.sh"
fi

{
  plumed config makefile_conf
  echo "PLUMED_INCLUDE=-I$(plumed info --include-dir)/plumed"
  echo "PLUMED_KERNEL=-L${PLUMED_KERNEL}"
} > Make.inc
