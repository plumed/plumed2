#!/usr/bin/env bash
# Ensure a kernel can be loaded by a host that supplies neither OpenMP nor MPI.
set -euo pipefail

prefix=${1:?Usage: check-kernel.sh PLUMED_PREFIX}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/plumed-kernel-load.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
compile_args=("$tmp/load.c" -o "$tmp/load")
if [[ $(uname) == Darwin ]]; then
  extension=dylib
else
  extension=so
  compile_args+=(-ldl)
fi

cat > "$tmp/load.c" <<'C'
#include <dlfcn.h>
#include <stdio.h>

int main(int argc, char **argv) {
  void *kernel;
  if(argc != 2) return 2;
  kernel = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if(!kernel) {
    fprintf(stderr, "Cannot load standalone PLUMED kernel: %s\n", dlerror());
    return 1;
  }
  return dlclose(kernel) != 0;
}
C
"${CC:-cc}" "${compile_args[@]}"
"$tmp/load" "$prefix/lib/libplumedKernel.$extension"
echo 'Standalone kernel loaded successfully from a plain C host'
