#!/usr/bin/env bash
# Build the selected serial or MPI sander executable.
set -euo pipefail

source_dir=${1:?Usage: build-amber.sh SOURCE BUILD}
build_dir=${2:?Usage: build-amber.sh SOURCE BUILD}
engine_mode=${MD_ENGINE_MODE:-mpi}
if [[ "$engine_mode" == mpi ]]; then
  : "${MD_MPICC:?Set MD_MPICC to the selected MPI C wrapper}"
  : "${MD_MPICXX:?Set MD_MPICXX to the selected MPI C++ wrapper}"
  : "${MD_MPIFC:?Set MD_MPIFC to the selected MPI Fortran wrapper}"
  mpi=TRUE
  target=sander.MPI
elif [[ "$engine_mode" == serial ]]; then
  mpi=FALSE
  target=sander
else
  echo "Unsupported MD_ENGINE_MODE: $engine_mode" >&2
  exit 2
fi

# MPI discovery must probe the same compilers that build Amber, including the
# Fortran module ABI. This also supports conda wrappers whose default driver
# names are unavailable outside their original build environment.
if [[ "$engine_mode" == mpi ]]; then
  export MPICH_CC="${CC:-cc}" MPICH_CXX="${CXX:-c++}" MPICH_FC="${FC:-gfortran}"
  export OMPI_CC="$MPICH_CC" OMPI_CXX="$MPICH_CXX" OMPI_FC="$MPICH_FC"
fi

cmake_args=(
  "-DCMAKE_C_COMPILER=${CC:-cc}"
  "-DCMAKE_CXX_COMPILER=${CXX:-c++}"
  "-DCMAKE_Fortran_COMPILER=${FC:-gfortran}"
  -DCOMPILER=MANUAL "-DMPI=$mpi" -DCUDA=FALSE
  -DBUILD_PYTHON=FALSE -DCHECK_UPDATES=FALSE -DINSTALL_TESTS=FALSE
  -DDOWNLOAD_MINICONDA=FALSE -DFORCE_DISABLE_LIBS=boost
  -DTRUST_SYSTEM_LIBS=TRUE -DCMAKE_BUILD_TYPE=Release
  "-DCMAKE_INSTALL_PREFIX=$build_dir/install"
)
if [[ "$engine_mode" == mpi ]]; then
  cmake_args+=("-DMPI_C_COMPILER=$MD_MPICC"
    "-DMPI_CXX_COMPILER=$MD_MPICXX" "-DMPI_Fortran_COMPILER=$MD_MPIFC")
fi

cmake -S "$source_dir" -B "$build_dir" \
  "${cmake_args[@]}"
cmake --build "$build_dir" --target "$target" --parallel "${MD_BUILD_JOBS:-4}"
