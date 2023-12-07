#! /usr/bin/env bash
#assign a default value to python bin and plumed
python_bin=${python_bin:-python}
plumed_program_name=${plumed_program_name:-plumed}
#python3-config may be any python3.* version, this should avoid contamination from other environments
#and use pythonX.Y-config, as suggested in the manual
pyver=$($python_bin -c "import sysconfig;print(sysconfig.get_config_var('VERSION'))")
python_config=python${pyver}-config

# Note that in conda libpython3xx is not found in the path returned by ldflags. IMHO it is a bug.
# The workaround is to -L appropriately. Will be fixed here.
conda_fixup=${CONDA_PREFIX+-L$CONDA_PREFIX/lib}
if [ ! -z "$conda_fixup" ]; then
  echo "CONDA_PREFIX is set. Assuming conda and enabling a workaround for missing -L in ${python_config} --ldflags --embed"
fi

if ! ${python_config} --embed >/dev/null 2>/dev/null; then
  #TODO: verify that this does not give problems with conda
  echo "PyCV needs python to be built to be embedable"
  echo "(compiling python with --enable-shared should be enough)"
  exit 1
fi

#for some reason, on the WSL I need to compile with `export PYCV_EXTRA_LDFLAGS="-Wl,--no-as-needed"`
#-Wl,--no-as-needed forces the python library to be linked, without this in a WSL does not work

#-fvisibility=hidden is needed to correct the warnings for the visibility of some pybind11 functionalities
export PLUMED_MKLIB_CFLAGS="$(${python_config} --cflags --embed) $(${python_bin} -m pybind11 --includes) -fvisibility=hidden"

export PLUMED_MKLIB_LDFLAGS="${PYCV_EXTRA_LDFLAGS} $(${python_config} --ldflags --embed) $conda_fixup"

echo PLUMED_MKLIB_CFLAGS=$PLUMED_MKLIB_CFLAGS
echo PLUMED_MKLIB_LDFLAGS=$PLUMED_MKLIB_LDFLAGS

${plumed_program_name} mklib PythonCVInterface.cpp ActionWithPython.cpp PythonFunction.cpp PlumedPythonEmbeddedModule.cpp
