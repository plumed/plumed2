#!/usr/bin/env bash
set -euo pipefail
export PLUMED_CONDA_BUILD=1
# The generic extension route (extensions/build.py), as for the other modules:
# the .cu files are compiled by NVCC with the flags of the base (src/maketools/cuda-flags),
# the library is linked by the C++ compiler with the shared CUDA runtime, and registered
# with the variant cuda<version of NVCC>.
nvcc=${CUDACXX:-nvcc}
actual_cuda=$("$nvcc" --version | sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\),.*/\1/p')
echo "Building CUDA variant $cuda_compiler_version with NVCC $actual_cuda"
[[ $actual_cuda == "$cuda_compiler_version" ]]
# Use the architecture matrix from conda-forge's compiler activation. For older
# activation packages, query NVCC rather than silently selecting one old GPU.
gencode=${NVCC_GENCODE:-}
if [[ -z $gencode ]]; then
  while read -r sm; do
    [[ $sm =~ ^sm_[0-9]+$ ]] || continue
    gencode+=" -gencode=arch=compute_${sm#sm_},code=$sm"
  done < <("$nvcc" --list-gpu-code)
  latest=$("$nvcc" --list-gpu-arch | awk '/^compute_[0-9]+$/' | sort -V | tail -n 1)
  gencode+=" -gencode=arch=$latest,code=$latest"
fi
NVCC="$nvcc" NVCCFLAGS="$gencode" ./configure --with-extension-base="$PREFIX" \
  --enable-extension-modules=gpu_cuda_hip --with-cudart=shared --prefix="$PREFIX"
make
make install
