#!/usr/bin/env bash
set -euo pipefail
ulimit -c 0
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
[[ $(uname) == Linux ]] || exit 0
# These checks inspect loading and package selection, not GPU calculations.
export CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA:-13.0}
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-cuda-install.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
cp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
for cuda_version in 12.9 13.0; do
  conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
    "plumed-ext-cuda-coordination=$VERSION" "cuda-version=$cuda_version"
  cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
  python extensions/tests/abi.py --prefix "$prefix"
  "$prefix/bin/plumed" modules -q --has cuda_coordination
  PLUMED_AUTOLOAD=0 "$prefix/bin/plumed" --load-module=cuda_coordination modules -q --has-feature cuda
  grep -qx "variant cuda$cuda_version" "$prefix/lib/plumed/plugins/cuda_coordination.plumed"
  for record in "$prefix"/conda-meta/{*nvcc*,cuda-compiler,python,py-plumed}-*.json; do
    test ! -e "$record"
  done
  if [[ $cuda_version == 12.9 ]]; then
    mkdir -p "$task/older/lib/plumed/plugins"
    cp "$prefix/lib/CudaCoordination.so" "$task/older/lib/"
    cp "$prefix/lib/plumed/plugins/cuda_coordination.plumed" "$task/older/lib/plumed/plugins/"
  else
    # A separately installed old variant must be rejected before native loading.
    if PLUMED_PLUGIN_PATH="$task/older/lib/plumed/plugins" \
      "$prefix/bin/plumed" modules -q --has cuda_coordination > "$task/ambiguous.log" 2>&1; then
      echo 'Loader accepted two CUDA coordination variants' >&2
      exit 1
    fi
    grep -q 'Ambiguous module cuda_coordination' "$task/ambiguous.log"
  fi
done
conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-ext-cuda-coordination
cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
! "$prefix/bin/plumed" modules -q --has cuda_coordination
"$prefix/bin/plumed" --has-thread-mpi
