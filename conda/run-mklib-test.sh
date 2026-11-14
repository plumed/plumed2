#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the built package version}"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-mklib-install.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
kernel=$(find "$prefix/lib" -maxdepth 1 -name 'libplumedKernel.*' -type f | head -n 1)
cp "$kernel" "$task/kernel-before"
test ! -f "$prefix/lib/plumed/scripts/mklib.sh"
if "$prefix/bin/plumed" mklib --help > "$task/missing.log" 2>&1; then
  echo 'Base package unexpectedly includes runtime compilation tools' >&2
  exit 1
fi
grep -q 'install plumed-mklib' "$task/missing.log"
conda install -y -p "$prefix" --override-channels -c local -c conda-forge "plumed-mklib=$VERSION"
cmp "$kernel" "$task/kernel-before"
conda install -y -p "$prefix" --override-channels -c local -c conda-forge "plumed-ext-ves=$VERSION" py-plumed python
if [[ $(uname) == Darwin ]]; then
  # conda hardlinks the toolchain executables to its package cache. dyld can then
  # resolve ld64's @loader_path/../lib through the cache hardlink and miss
  # libtapi.dylib (which lives in this environment, not in ld64's cache folder).
  # Copy the exact toolchain packages selected above, as for openmpi in
  # run-conda-test.sh.
  toolchain=$(conda list --prefix "$prefix" --export | grep -E '^(ld64|cctools|clang)')
  conda install -y --copy --force-reinstall --freeze-installed -p "$prefix" \
    --override-channels -c local -c conda-forge $toolchain
fi
conda run --no-capture-output -p "$prefix" python extensions/tests/mklib.py --prefix "$prefix" --extension
"$prefix/bin/plumed" mklib --help
conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-mklib
cmp "$kernel" "$task/kernel-before"
test ! -f "$prefix/lib/plumed/scripts/mklib.sh"
test ! -f "$prefix/lib/plumed/plumed-mklib"
