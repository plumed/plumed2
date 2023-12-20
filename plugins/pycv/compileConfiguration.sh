#!/bin/env bash

#formatted with shfmt  https://github.com/mvdan/sh/releases
#checked with shellcheck
# shellcheck disable=SC2034  # The unused variables will be used in sourcer files

#please source this file

#for some reason, on the WSL I need to compile with 
# `export PYCV_EXTRA_LDFLAGS="-Wl,--no-as-needed"`
#-Wl,--no-as-needed forces the python library to be linked
#I do not undestand why, since -as-needed should be disabled by default

plumed_program_name=${plumed_program_name:-plumed}

#check if during config time the settings for compiling PyCV where got
plumed_canPyCV=$(
  ${plumed_program_name} --no-mpi config makefile_conf | grep canPyCV | sed -e s/^canPyCV=//
)

if [[ $plumed_canPyCV = yes ]]; then
  pybind11_cflags=$(plumed --no-mpi config makefile_conf | grep pybind11_cflags | sed -e s/^pybind11_cflags=//)
  python_cf_embedded=$(plumed --no-mpi config makefile_conf | grep python_cf_embedded | sed -e s/^python_cf_embedded=//)
  python_ld_embedded=$(plumed --no-mpi config makefile_conf | grep python_ld_embedded | sed -e s/^python_ld_embedded=//)
else
  #assign a default value to python bin and plumed
  python_bin=${python_bin:-python}
  #python3-config may be any python3.* version, this should avoid contamination from other environments
  #and use pythonX.Y-config, as suggested in the manual
  pyver=$($python_bin -c "import sysconfig;print(sysconfig.get_config_var('VERSION'))")
  python_config=python${pyver}-config

  # Note that in conda libpython3xx is not found in the path returned by ldflags. IMHO it is a bug.
  # The workaround is to -L appropriately. Will be fixed here.
  conda_fixup=${CONDA_PREFIX+-L$CONDA_PREFIX/lib}
  if [ -n "$conda_fixup" ]; then
    echo "CONDA_PREFIX is set. Assuming conda and enabling a workaround for missing -L in ${python_config} --ldflags --embed"
  fi

  if ! ${python_config} --embed >/dev/null 2>/dev/null; then
    #TODO: verify that this does not give problems with conda
    echo "PyCV needs python to be built to be embedable"
    echo "(compiling python with --enable-shared should be enough)"
    exit 1
  fi
  if ! ${python_bin} -m pybind11 --includes >/dev/null; then
    #python will automatically say that there is non pybind11
    exit 1
  fi
  #-fvisibility=hidden is needed to correct the warnings for the visibility of some pybind11 functionalities
  pybind11_cflags="$(${python_bin} -m pybind11 --includes) -fvisibility=hidden"
  python_cf_embedded=$(${python_config} --cflags --embed)
  python_ld_embedded=$(${python_config} --ldflags --embed)
fi
