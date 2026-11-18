#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-pycv-install.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
kernel=$(find "$prefix/lib" -maxdepth 1 -name 'libplumedKernel.*' -type f | head -n 1)
cp "$kernel" "$task/kernel-before"
for python_version in 3.11 3.15; do
  conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
    "plumed-ext-pycv=$VERSION" "python=$python_version"
  cmp "$kernel" "$task/kernel-before"
  "$prefix/bin/python" -c "import sys; assert '.'.join(map(str, sys.version_info[:2])) == '$python_version'"
  "$prefix/bin/plumed" modules -q --has pycv
  PLUMED_AUTOLOAD=0 "$prefix/bin/plumed" --load-module=pycv modules -q --has pycv
  # Only the plugin needs Python; neither wrappers nor build tools are required.
  for record in "$prefix"/conda-meta/{py-plumed,pybind11,cmake}-*.json; do
    test ! -e "$record"
  done
  conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-ext-pycv
  cmp "$kernel" "$task/kernel-before"
  ! "$prefix/bin/plumed" modules -q --has pycv
  "$prefix/bin/plumed" --has-thread-mpi
done
