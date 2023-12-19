#! /usr/bin/env bash

#formatted with shfmt  https://github.com/mvdan/sh/releases
#checked with shellcheck
# the SC2154 warnings are variables present in the sourced file

source compileConfiguration.sh


if [[ -z $PLUMED_KERNEL ]]; then
  echo "$(basename $0) can work only if \"PLUMED_KERNEL\" is defined"
  echo "either via module load or sourceme.sh, or manually exported"
  exit 1
fi

plumed_include=$(${plumed_program_name} --no-mpi info --include-dir)
{
  ${plumed_program_name} --no-mpi config makefile_conf
  echo "PLUMED_KERNEL=${PLUMED_KERNEL}"
  if [[ ! $plumed_canPyCV = yes ]]; then
    echo "pybind11_cflags=$pybind11_cflags"
    echo "python_cf_embedded=$python_cf_embedded"
    echo "python_ld_embedded=$python_ld_embedded"
  fi
  echo "ADDCPPFLAGS=-I${plumed_include}/plumed"
  echo "ADDCLDFLAGS=${PYCV_EXTRA_LDFLAGS} ${conda_fixup}"
} > Makefile.conf
