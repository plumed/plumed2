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

ls -l $CONDA_PREFIX/conda-bld/
ls -l $CONDA_PREFIX/conda-bld/*
