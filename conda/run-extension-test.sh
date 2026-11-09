#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-extension-test.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
kernel=$(find "$prefix/lib" -maxdepth 1 -name 'libplumedKernel.*' -type f | head -n 1)
cp "$kernel" "$task/kernel-before"
! "$prefix/bin/plumed" config -q module ves
conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
  "plumed-ext-ves=$VERSION" "plumed-ext-pamm=$VERSION"
cmp "$kernel" "$task/kernel-before"
python extensions/tests/installed.py --base "$prefix" --extensions "$prefix"
conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-ext-ves
cmp "$kernel" "$task/kernel-before"
! "$prefix/bin/plumed" config -q module ves
"$prefix/bin/plumed" config -q module pamm adjmat
# Verify removal leaves all base functionality intact.
"$prefix/bin/plumed" --has-thread-mpi
"$prefix/bin/plumed" config -q module bias
