#!/usr/bin/env bash
set -euo pipefail
export PLUMED_AUTOLOAD=0
export PLUMED_KERNEL="$PREFIX/lib/libplumedKernel.so"
python extensions/register.py --plumed "$PREFIX/bin/plumed" --module cuda_coordination \
  --library ../../CudaCoordination.so --variant "cuda${cuda_compiler_version}" \
  --features CUDA --output "$SRC_DIR/cuda-registration"
"$CXX" $CXXFLAGS $CPPFLAGS -std=c++17 -fPIC -D__PLUMED_HAS_CUDA=1 \
  -I"$PREFIX/include/plumed" -c "$SRC_DIR/cuda-registration/module_identity.cpp" \
  -o "$SRC_DIR/cuda-registration/module_identity.o"
cd plugins/cudaCoord
bash configure.sh
# Use the architecture matrix from conda-forge's compiler activation. For older
# activation packages, query NVCC rather than silently selecting one old GPU.
nvcc=${CUDACXX:-nvcc}
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
make -j"${CPU_COUNT:-4}" CXX="$CXX" NVCC="$nvcc" \
  EXTRA_OBJS="$SRC_DIR/cuda-registration/module_identity.o" \
  NVCCCFLAGS="-dc --forward-unknown-to-host-compiler $gencode" \
  NVCCLDFLAGS="-shared --cudart=shared --forward-unknown-to-host-compiler $gencode"
install -m 755 CudaCoordination.so "$PREFIX/lib/"
mkdir -p "$PREFIX/lib/plumed/plugins"
cp "$SRC_DIR/cuda-registration/cuda_coordination.plumed" "$PREFIX/lib/plumed/plugins/"
