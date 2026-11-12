#! /usr/bin/env bash

set -e
cd "$(dirname "$0")"

# macOS dyld can resolve hardlinked compiler paths through the package cache.
# Copy toolchains into each build prefix so their relative dylibs stay together.
if [[ $(uname) == Darwin ]]; then
  export CONDA_ALWAYS_COPY=true
fi

export CPU_COUNT=${CPU_COUNT:-4}
# Resolve unpinned Python-wrapper tests against this CI build first.
export CONDA_CHANNEL_PRIORITY=strict
echo "#Building plumed for conda"
conda-build -c conda-forge plumed
# Prefer the package just built locally, not an existing conda-forge PLUMED.
echo "#Building optional action modules"
python ../extensions/conda_recipe.py
conda-build -c local -c conda-forge plumed-extensions
echo "#Building the MPICH and Open MPI add-ons"
conda-build -c local -c conda-forge plumed-mpi
echo "#Building py-plumed for conda"
conda-build -c local -c conda-forge py-plumed
echo "#Building runtime C++ compilation tools"
conda-build -c local -c conda-forge plumed-mklib
echo "#Building Python collective-variable plugins"
conda-build -c local -c conda-forge plumed-pycv
if [[ $(uname) == Linux ]]; then
  echo "#Building CUDA coordination variants"
  # CI has no GPU: override driver metadata for compilation/loading tests only.
  CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA:-13.0} \
    conda-build -c local -c conda-forge plumed-cuda-coordination
fi

if [[ $(uname) == Linux && $(uname -m) == x86_64 && -n ${NVHPC_HOME:-} ]]; then
  echo "#Building the multi-target OpenACC plugin"
  conda-build -c local -c conda-forge plumed-openacc
fi

build_root=${CONDA_BLD_PATH:-$CONDA_PREFIX/conda-bld}
ls -l "$build_root"
ls -l "$build_root"/*
