#!/usr/bin/env bash
set -euo pipefail
ulimit -c 0
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
[[ $(uname) == Linux ]] || exit 0
# These checks inspect loading and package selection, not GPU calculations.
# gpu_voronoi uses only the headers of gpu_cuda_hip: its package must not pull
# plumed-ext-gpu-cuda-hip in, and works without it.
export CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA:-$(python conda/gpu-plan.py --max-cuda)}
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-gpu-voronoi-install.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
cp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
for cuda_version in $(python conda/gpu-plan.py --cuda-versions); do
  conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
    "plumed-ext-gpu-voronoi=$VERSION" "cuda-version=$cuda_version"
  cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
  for record in "$prefix"/conda-meta/plumed-ext-gpu-cuda-hip-*.json; do
    test ! -e "$record"
  done
  "$prefix/bin/plumed" modules -q --has gpu_voronoi
  if "$prefix/bin/plumed" modules -q --has gpu_cuda_hip ; then
    echo 'gpu_cuda_hip is installed with gpu_voronoi' >&2
    exit 1
  fi
  PLUMED_AUTOLOAD=0 "$prefix/bin/plumed" --load-module=gpu_voronoi modules -q --has-feature cuda
  grep -qx "variant cuda$cuda_version" "$prefix/lib/plumed/plugins/gpu_voronoi.plumed"
done
conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-ext-gpu-voronoi
cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
if "$prefix/bin/plumed" modules -q --has gpu_voronoi ; then
  echo 'gpu_voronoi is still installed' >&2
  exit 1
fi
