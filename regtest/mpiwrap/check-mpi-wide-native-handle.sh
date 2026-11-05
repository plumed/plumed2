#!/usr/bin/env bash
# Native MPI remains usable even when its handles exceed wrapper capacity.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-wide-native.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
cat > "$TASK_DIR/mpi.h" <<'HEADER'
#ifndef MPIW_TEST_MPI_H
#define MPIW_TEST_MPI_H
typedef struct { unsigned char bytes[16]; } MPI_Comm;
typedef void *MPI_Datatype;
typedef void *MPI_Errhandler;
typedef void *MPI_Group;
typedef void *MPI_Op;
typedef void *MPI_Request;
#endif
HEADER
cat > "$TASK_DIR/probe.cpp" <<'SOURCE'
#define __PLUMED_HAS_MPI 1
#include "mpiwrap.h"
#include <cassert>
#include <cstring>
namespace PLMD {
int checkWideNativeHandle() {
  MPI_Comm original;
  std::memset(&original, 0xa5, sizeof(original));
  const MPI_Comm imported = mpiwrap_Comm_from_native(&original);
  assert(std::memcmp(&original, &imported, sizeof(original)) == 0);
  unsigned char output[sizeof(MPI_Comm)];
  std::memset(output, 0x5a, sizeof(output));
  assert(!mpiwrap_Comm_export_native(imported, output, sizeof(output) - 1));
  for (unsigned char byte : output) assert(byte == 0x5a);
  assert(mpiwrap_Comm_export_native(imported, output, sizeof(output)));
  assert(std::memcmp(&original, output, sizeof(original)) == 0);
  return 0;
}
}
int main() { return PLMD::checkWideNativeHandle(); }
SOURCE
for capacity in 8 16; do
  "${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -pedantic \
    -DPLUMED_MPIW_HANDLE_CAPACITY=$capacity \
    -I"$TASK_DIR" -I"$ROOT/src/mpiwrap" "$TASK_DIR/probe.cpp" -o "$TASK_DIR/probe-$capacity"
  "$TASK_DIR/probe-$capacity"
done
echo '16-byte native communicator works independently of wrapper capacity (8 and 16)'
