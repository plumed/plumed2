#!/usr/bin/env bash
# Build an already patched GROMACS tree. Also usable outside GitHub Actions.
set -euo pipefail

source_dir=${1:?Usage: build-gromacs.sh SOURCE BUILD}
build_dir=${2:?Usage: build-gromacs.sh SOURCE BUILD}
: "${BACKEND:?Set BACKEND to tmpi, openmpi, or mpich}"
options=()
case "$BACKEND" in
  tmpi) options+=(-DGMX_MPI=OFF -DGMX_THREAD_MPI=ON) ;;
  openmpi|mpich)
    options+=(-DGMX_MPI=ON -DGMX_THREAD_MPI=OFF
      "-DMPI_C_COMPILER=${MD_MPICC:?Set MD_MPICC}"
      "-DMPI_CXX_COMPILER=${MD_MPICXX:?Set MD_MPICXX}")
    ;;
  *) echo "Unknown backend: $BACKEND" >&2; exit 2 ;;
esac

if [[ -n "${MD_LIBOMP:-}" ]]; then
  options+=("-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I$MD_LIBOMP/include"
    -DOpenMP_CXX_LIB_NAMES=omp "-DOpenMP_omp_LIBRARY=$MD_LIBOMP/lib/libomp.dylib")
fi

cmake -S "$source_dir" -B "$build_dir" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_C_COMPILER="${CC:-cc}" \
  -DCMAKE_CXX_COMPILER="${CXX:-c++}" -DGMX_GPU=OFF -DGMX_OPENMP=ON \
  -DGMX_USE_PLUMED=ON -DGMX_FFT_LIBRARY=fftpack -DGMXAPI=OFF \
  -DBUILD_TESTING=OFF -DGMX_BUILD_UNITTESTS=OFF -DGMX_BUILD_HELP=OFF \
  -DREGRESSIONTEST_DOWNLOAD=OFF "${options[@]}"
cmake --build "$build_dir" -j "${MD_BUILD_JOBS:-4}"
if [[ "$BACKEND" == tmpi ]]; then
  "$build_dir/bin/gmx" --version
else
  "$build_dir/bin/gmx_mpi" --version
fi
