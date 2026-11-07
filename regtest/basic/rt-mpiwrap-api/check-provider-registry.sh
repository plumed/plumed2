#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/../../.." && pwd)
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/mpiw-registry.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM

for capacity in default 16 24; do
  capacity_flag=
  if [ "$capacity" = 16 ]; then capacity_flag=-DPLUMED_MPIW_HANDLE_CAPACITY=16; fi
  if [ "$capacity" = 24 ]; then
    capacity_flag="-DPLUMED_MPIW_HANDLE_CAPACITY=24 -DPLUMED_MPIW_STATUS_CAPACITY=192"
  fi
  "${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -pthread -D__PLUMED_MPIWRAP $capacity_flag \
    -I"$repo_root/src" -iquote "$repo_root/src/mpiwrap" \
    -x c++ "$repo_root/src/mpiwrap/mpiwrap.cpp" \
    -x c++ "$repo_root/src/mpiwrap/loader.cpp" \
    -x c++ "$repo_root/src/mpiwrap/fallbacks.cpp" \
    -x c++ "$repo_root/src/mpiwrap/launcher.cpp" \
    -x c++ "$repo_root/src/mpiwrap/tmpi.cpp" \
    -x c++ "$repo_root/regtest/basic/rt-mpiwrap-api/test-provider-registry.cpp.in" \
    -o "$tmpdir/test-provider-registry-$capacity"
  "$tmpdir/test-provider-registry-$capacity"
done
echo "provider registry validation and typed restoration passed with default, 16-byte, and configured storage"
