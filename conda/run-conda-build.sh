#! /usr/bin/env bash

set -e
cd "$(dirname "$0")"

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
echo "#Building Python collective-variable plugins"
conda-build -c local -c conda-forge plumed-pycv
if [[ $(uname) == Linux ]]; then
  echo "#Building CUDA coordination variants"
  # CI has no GPU: override driver metadata for compilation/loading tests only.
  CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA:-13.0} \
    conda-build -c local -c conda-forge plumed-cuda-coordination
fi

ls -l $CONDA_PREFIX/conda-bld/
ls -l $CONDA_PREFIX/conda-bld/*
