#!/usr/bin/env bash
# Check the standalone wrapper without inheriting a configured kernel's flags.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
CXX=${CXX:-c++}
CC=${CC:-cc}

"$CC" -std=c99 -D__PLUMED_HAS_DLOPEN -fsyntax-only "$ROOT/src/wrapper/Plumed.c"
for standard in c++98 c++11 c++17; do
  for external in 0 1; do
    "$CXX" -std="$standard" -D__PLUMED_HAS_DLOPEN \
      -D__PLUMED_WRAPPER_IMPLEMENTATION=1 -D__PLUMED_WRAPPER_EXTERN="$external" \
      -x c++ -fsyntax-only "$ROOT/src/wrapper/Plumed.h"
  done
done
echo "standalone wrapper C99 and C++98/11/17 compatibility passed"
