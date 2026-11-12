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
if [[ ${PLUMED_CONDA_SKIP_BASE:-no} == yes ]]; then
  : "${CONDA_BLD_PATH:?set CONDA_BLD_PATH for the shared kernel archive}"
  : "${VERSION:?set VERSION to the shared kernel package version}"
  subdir=$(conda info --json | python -c 'import json,sys; print(json.load(sys.stdin)["platform"])')
  if ! compgen -G "$CONDA_BLD_PATH/$subdir/plumed-$VERSION-*.conda" >/dev/null; then
    echo "Missing shared kernel archive for $subdir/$VERSION" >&2
    exit 1
  fi
else
  echo "#Building plumed for conda"
  conda-build -c conda-forge plumed
fi
# Prefer the exact package built or restored locally to conda-forge PLUMED.
python ../extensions/conda_recipe.py
# Exercise runtime compilation before the more expensive module/plugin builds.
echo "#Building py-plumed for conda"
conda-build -c local -c conda-forge py-plumed
echo "#Building runtime C++ compilation tools"
conda-build -c local -c conda-forge plumed-mklib
echo "#Building optional action modules"
conda-build -c local -c conda-forge plumed-extensions
echo "#Building the MPICH and Open MPI add-ons"
conda-build -c local -c conda-forge plumed-mpi
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
