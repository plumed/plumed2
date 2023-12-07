#! /usr/bin/env bash

#assign a default value to python bin and plumed
python_bin=${python_bin:-python}
plumed_program_name=${plumed_program_name:-plumed}
#python3-config may be any python3.* version, this should avoid contamination from other environments
#and use pythonX.Y-config, as suggested in the manual
pyver=$($python_bin -c "import sysconfig;print(sysconfig.get_config_var('VERSION'))")
python_config=python${pyver}-config

if [[ -z $PLUMED_KERNEL ]]; then
  echo "$(basename $0) can work only if \"PLUMED_KERNEL\" is defined"
  echo "either via module load or sourceme.sh, or manually exported"
  exit 1
fi

if ! ${python_config} --embed >/dev/null 2>/dev/null; then
#TODO: verify that this does not give problems with conda
  echo "PyCV needs python to be built to be embedable"
  echo "(compiling python with --enable-shared should be enough)"
  exit 1
fi

{
  ${plumed_program_name} --no-mpi config makefile_conf
  echo "PLUMED_KERNEL=${PLUMED_KERNEL}"
  echo "ADDCPPFLAGS=$(${python_config} --cflags --embed) $(${python_bin} -m pybind11 --includes) -I$(${plumed_program_name} --no-mpi info --include-dir)/plumed"
  echo "ADDCLDFLAGS=$(PYCV_EXTRA_LDFLAGS) $(${python_config} --ldflags --embed)"
} > Make.inc
