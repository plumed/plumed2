#! /usr/bin/env bash

if [[ -z $PLUMED_KERNEL ]]; then
  echo "$(basename $0) can work only if \"PLUMED_KERNEL\" is defined"
  echo "either via module load or sourceme.sh, or manually exported"
  exit 1
fi

if ! python3-config --embed >/dev/null 2>/dev/null; then
#TODO: verify that this does not give problems with conda
  echo "PyCV needs python to be built to be embedable"
  echo "(compiling python with --enable-shared should be enough)"
  exit 1
fi

{
  plumed --no-mpi config makefile_conf
  echo "PLUMED_KERNEL=-L${PLUMED_KERNEL}"
  echo "ADDCPPFLAGS=$(python3-config --cflags --embed) $(python3 -m pybind11 --includes) -I$(plumed --no-mpi info --include-dir)/plumed"
  echo "ADDCLDFLAGS=$(python3-config --ldflags --embed)"
} > Make.inc
