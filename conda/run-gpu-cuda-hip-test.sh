#!/usr/bin/env bash
set -euo pipefail
ulimit -c 0
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
[[ $(uname) == Linux ]] || exit 0
# These checks inspect loading and package selection, not GPU calculations.
export CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA:-$(python conda/gpu-plan.py --max-cuda)}
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-gpu-cuda-hip-install.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
cp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
first_variant=yes
for cuda_version in $(python conda/gpu-plan.py --cuda-versions); do
  conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
    "plumed-ext-gpu-cuda-hip=$VERSION" "cuda-version=$cuda_version"
  cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
  python extensions/tests/abi.py --prefix "$prefix"
  "$prefix/bin/plumed" modules -q --has gpu_cuda_hip
  PLUMED_AUTOLOAD=0 "$prefix/bin/plumed" --load-module=gpu_cuda_hip modules -q --has-feature cuda
  grep -qx "variant cuda$cuda_version" "$prefix/lib/plumed/plugins/gpu_cuda_hip.plumed"
  for record in "$prefix"/conda-meta/{*nvcc*,cuda-compiler,python,py-plumed}-*.json; do
    test ! -e "$record"
  done
  if [[ $first_variant == yes ]]; then
    first_variant=no
    mkdir -p "$task/older/lib/plumed/plugins"
    cp "$prefix/lib/libplumed-ext-gpu_cuda_hip.so" "$task/older/lib/"
    cp "$prefix/lib/plumed/plugins/gpu_cuda_hip.plumed" "$task/older/lib/plumed/plugins/"
  else
    # A separately installed old variant must be rejected before native loading.
    if PLUMED_PLUGIN_PATH="$task/older/lib/plumed/plugins" \
      "$prefix/bin/plumed" modules -q --has gpu_cuda_hip > "$task/ambiguous.log" 2>&1; then
      echo 'Loader accepted two gpu_cuda_hip variants' >&2
      exit 1
    fi
    grep -q 'Ambiguous module gpu_cuda_hip' "$task/ambiguous.log"
  fi
done
conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-ext-gpu-cuda-hip
cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
! "$prefix/bin/plumed" modules -q --has gpu_cuda_hip
"$prefix/bin/plumed" --has-thread-mpi
