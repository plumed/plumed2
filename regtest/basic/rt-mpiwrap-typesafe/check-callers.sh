#!/usr/bin/env bash
# C11 callers need no MPI runtime: exercise both legal native representations.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
fixture="$root/regtest/basic/rt-mpiwrap-typesafe/caller.c.in"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-mpi-typesafe.XXXXXX")
trap 'rm -rf "$task"' EXIT
for representation in int 'struct opaque_comm *' 'void *'; do
  flag=
  if [[ $representation = int ]]; then flag=-DPLUMED_TEST_INT_COMM; fi
  if [[ $representation = 'void *' ]]; then flag=-DPLUMED_TEST_VOID_COMM; fi
  cat > "$task/mpi.h" <<HEADER
#ifndef MPI_TYPESAFE_TEST_H
#define MPI_TYPESAFE_TEST_H
#define MPI_VERSION 5
typedef $representation MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Errhandler;
typedef int MPI_Group;
typedef int MPI_Op;
typedef int MPI_Request;
#endif
HEADER
  for ccheck in 0 1; do
    "${CC:-cc}" -D__PLUMED_WRAPPER_C_MPI_COMM_CHECK=$ccheck -std=c11 -Wall -Wextra -Werror $flag -I"$task" -x c "$fixture" -o "$task/caller"
    "$task/caller"
  done
  "${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -Wno-unused-parameter -D__PLUMED_HAS_MPI \
    -I"$task" -I"$root/src" -x c++ "$root/regtest/basic/rt-mpiwrap-typesafe/kernel.cpp.in" \
    "$root/src/tools/TypesafePtr.cpp" "$root/src/tools/Exception.cpp" -o "$task/kernel"
  "$task/kernel"
  for standard in c++98 c++11 c++17; do
    for options in "" "-D__PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE=1" "-D__PLUMED_WRAPPER_CXX_MPI_COMM_CHECK=1" "-D__PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE=1 -D__PLUMED_WRAPPER_CXX_MPI_COMM_CHECK=0"; do
      "${CXX:-c++}" $options -std="$standard" -Wall -Wextra -Werror $flag -I"$task" -x c++ \
        "$root/regtest/basic/rt-mpiwrap-typesafe/caller.cpp.in" -o "$task/caller-cxx"
      "$task/caller-cxx"
    done
  done
done
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -DPLUMED_TEST_WRAPPER_COMM \
  -x c++ "$root/regtest/basic/rt-mpiwrap-typesafe/caller.cpp.in" -o "$task/caller-wrapper"
"$task/caller-wrapper"
# Type-only consumers must retain the same communicator annotations.
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -DPLUMED_TEST_WRAPPER_COMM -DPLUMED_TEST_TYPES_ONLY \
  -x c++ "$root/regtest/basic/rt-mpiwrap-typesafe/caller.cpp.in" -o "$task/caller-types"
"$task/caller-types"
# These common headers must not pull provider dispatch into every action.
cat > "$task/common.cpp" <<HEADER
#include "$root/src/tools/TypesafePtr.h"
#include "$root/src/tools/Communicator.h"
#ifdef __PLUMED_mpiwrap_mpiwrap_h
#error "common type consumers must not include the package interface"
#endif
int main() { return sizeof(MPI_Comm) == 0; }
HEADER
# Match normal non-OpenACC builds: the included headers contain acc pragmas.
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -Wno-unused-parameter -Wno-unknown-pragmas \
  -D__PLUMED_MPIWRAP -I"$root/src" -fsyntax-only "$task/common.cpp"

bash "$root/regtest/basic/rt-mpiwrap-typesafe/check-header-odr.sh"
