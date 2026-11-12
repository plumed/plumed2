#!/usr/bin/env bash
set -euo pipefail
ulimit -c 0
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-openacc-install.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
# No driver metadata override: the multi-target package must install without GPU.
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
cp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
conda install -y -p "$prefix" --override-channels -c local -c conda-forge "plumed-ext-openacc-ptm=$VERSION"
cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
python extensions/tests/abi.py --prefix "$prefix"
NVCOMPILER_ACC_DEVICE_TYPE=HOST "$prefix/bin/plumed" modules -q --has openacc_ptm
PLUMED_AUTOLOAD=0 NVCOMPILER_ACC_DEVICE_TYPE=HOST \
  "$prefix/bin/plumed" --load-module=openacc_ptm modules -q --has-feature openacc
for record in "$prefix"/conda-meta/{*nvcc*,cuda-compiler,python,py-plumed,nvhpc}-*.json; do
  test ! -e "$record"
done
conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-ext-openacc-ptm
cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
! "$prefix/bin/plumed" modules -q --has openacc_ptm
"$prefix/bin/plumed" --has-thread-mpi
