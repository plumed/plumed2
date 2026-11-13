#!/usr/bin/env bash
set -euo pipefail
export PLUMED_AUTOLOAD=0
export PLUMED_KERNEL="$PREFIX/lib/libplumedKernel.so"
python extensions/register.py --plumed "$PREFIX/bin/plumed" --module cuda_coordination \
  --library ../../CudaCoordination.so --variant "cuda${cuda_compiler_version}" \
  --features CUDA --output "$SRC_DIR/cuda-registration"
# configure.sh translates the kernel flags for NVCC. Use the plugin's explicit
# CUDA rule; make's implicit C++ rule would pass -Xcompiler flags to plain g++.
cp "$SRC_DIR/cuda-registration/module_identity.cpp" "$SRC_DIR/plugins/cudaCoord/src/ModuleIdentity.cu"
policy="$SRC_DIR/extensions/conda_config.py"
cd plugins/cudaCoord
python "$policy" --module cuda_coordination -- bash configure.sh
# Use the architecture matrix from conda-forge's compiler activation. For older
# activation packages, query NVCC rather than silently selecting one old GPU.
nvcc=${CUDACXX:-nvcc}
actual_cuda=$("$nvcc" --version | sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\),.*/\1/p')
echo "Building CUDA variant $cuda_compiler_version with NVCC $actual_cuda"
[[ $actual_cuda == "$cuda_compiler_version" ]]
gencode=${NVCC_GENCODE:-}
if [[ -z $gencode ]]; then
  while read -r sm; do
    [[ $sm =~ ^sm_[0-9]+$ ]] || continue
    gencode+=" -gencode=arch=compute_${sm#sm_},code=$sm"
  done < <("$nvcc" --list-gpu-code)
  latest=$("$nvcc" --list-gpu-arch | awk '/^compute_[0-9]+$/' | sort -V | tail -n 1)
  gencode+=" -gencode=arch=$latest,code=$latest"
fi
# Fat binaries include all configured GPU targets. Device LTO is left off so
# each target carries normal device code through the separate link step.
python "$policy" --module cuda_coordination -- make -j"${CPU_COUNT:-4}" CXX="$CXX" NVCC="$nvcc" \
  PLUMED_INCLUDE="-I$PREFIX/include -I$PREFIX/include/plumed" \
  EXTRA_OBJS=ModuleIdentity.o \
  DYNAMIC_LIBS="-L$PREFIX/lib" \
  NVCCCFLAGS="-dc --forward-unknown-to-host-compiler -Xcompiler -pthread $gencode" \
  NVCCLDFLAGS="-shared --cudart=shared --forward-unknown-to-host-compiler -Xcompiler -pthread -Xlinker --no-undefined $gencode"
install -m 755 CudaCoordination.so "$PREFIX/lib/"
mkdir -p "$PREFIX/lib/plumed/plugins"
cp "$SRC_DIR/cuda-registration/cuda_coordination.plumed" "$PREFIX/lib/plumed/plugins/"

python "$policy" --module cuda_coordination --check-features --prefix "$PREFIX"
