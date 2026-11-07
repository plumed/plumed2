#!/usr/bin/env bash
# Build an already patched LAMMPS tree against the serial PLUMED installation.
set -euo pipefail

source_dir=${1:?Usage: build-lammps.sh SOURCE BUILD PLUMED_PREFIX}
build_dir=${2:?Usage: build-lammps.sh SOURCE BUILD PLUMED_PREFIX}
plumed_prefix=${3:?Usage: build-lammps.sh SOURCE BUILD PLUMED_PREFIX}
engine_mode=${MD_ENGINE_MODE:-mpi}
case "$engine_mode" in mpi|serial) ;; *) exit 2 ;; esac
build_mpi=ON
if [[ "$engine_mode" == serial ]]; then
  build_mpi=OFF
fi
cmake_args=(-DCMAKE_BUILD_TYPE=Release "-DCMAKE_CXX_COMPILER=${CXX:-c++}" "-DBUILD_MPI=$build_mpi")
if [[ "$engine_mode" == mpi ]]; then
  cmake_args+=("-DMPI_CXX_COMPILER=${MD_MPICXX:?Set MD_MPICXX}")
fi
if [[ -n "${MD_LIBOMP:-}" ]]; then
  cmake_args+=("-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I$MD_LIBOMP/include"
    -DOpenMP_CXX_LIB_NAMES=omp "-DOpenMP_omp_LIBRARY=$MD_LIBOMP/lib/libomp.dylib")
fi

export PKG_CONFIG_PATH="$plumed_prefix/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
cmake -S "$source_dir/cmake" -B "$build_dir" -G Ninja \
  "${cmake_args[@]}" \
  -DBUILD_OMP=ON -DPKG_OPENMP=ON -DPKG_PLUMED=ON \
  -DPLUMED_MODE=runtime -DDOWNLOAD_PLUMED=OFF -DDOWNLOAD_POTENTIALS=OFF \
  -DBUILD_SHARED_LIBS=OFF
cmake --build "$build_dir" -j "${MD_BUILD_JOBS:-4}"
"$build_dir/lmp" -help
