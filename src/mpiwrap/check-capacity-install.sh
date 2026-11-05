#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: check-capacity-install.sh MPICH_MPICXX OPENMPI_MPICXX" >&2
  exit 2
fi
mpich_mpicxx=$1
openmpi_mpicxx=$2
script_dir=$(cd "$(dirname "$0")" && pwd)
repo_root=$(cd "$script_dir/../.." && pwd)
if [[ ! -r "$repo_root/Makefile.conf" ]]; then
  echo "check-capacity-install.sh: configure the source tree first" >&2
  exit 2
fi
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiw-capacity-install.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

# Installed headers carry defaults; consumers choose matching overrides with -D.
capacity_flags=(-DPLUMED_MPIW_HANDLE_CAPACITY=24 -DPLUMED_MPIW_STATUS_CAPACITY=192)
make -C "$repo_root/src/lib" install-headers \
  DESTDIR="$tmp/stage" includedir=/usr/include program_name=plumed \
  DIRS=../mpiwrap

include_dir="$tmp/stage/usr/include"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror "${capacity_flags[@]}" \
  -I"$include_dir" -fsyntax-only -x c++ \
  "$script_dir/test-capacity-frontend.cpp.in"

MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 "${capacity_flags[@]}" \
  -D__PLUMED_HAS_MPI -DGMX_THREAD_MPI=0 \
  -I"$include_dir" -fsyntax-only -x c++ \
  "$script_dir/test-capacity-native.cpp.in"
OMPI_CXX="${OMPI_CXX:-c++}" "$openmpi_mpicxx" -std=c++17 "${capacity_flags[@]}" \
  -D__PLUMED_HAS_MPI -DGMX_THREAD_MPI=0 \
  -I"$include_dir" -fsyntax-only -x c++ \
  "$script_dir/test-capacity-native.cpp.in"

echo "compiler overrides give installed public headers 24-byte handles and 192-byte status storage in wrapper and native-MPI builds"
