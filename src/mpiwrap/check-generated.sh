#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
  echo "usage: check-generated.sh MPICH_MPICXX OPENMPI_MPICXX THREAD_MPI_INCLUDE" >&2
  exit 2
fi
mpich_mpicxx=$1
openmpi_mpicxx=$2
tmpi_include=$3
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/../.." && pwd)
generated="$repo_root/src/mpiwrap"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiw-generated.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

python3 "$repo_root/src/mpiwrap/generate.py" --check

"${CC:-cc}" -E -I"$generated" -DGMX_THREAD_MPI=1 -DGMX_VERSION=20260000 \
  -x c "$generated/test-gromacs-version-guard.c.in" >/dev/null
if "${CC:-cc}" -E -I"$generated" -DGMX_THREAD_MPI=1 \
    -x c "$generated/test-gromacs-version-guard.c.in" >"$tmp/no-gromacs-version.log" 2>&1; then
  echo "GROMACS tMPI adapter unexpectedly accepted an unknown host version" >&2
  exit 1
fi
grep -q 'requires GMX_VERSION or PLUMED_MPIW_GROMACS_VERSION' "$tmp/no-gromacs-version.log"
if "${CC:-cc}" -E -I"$generated" -DGMX_THREAD_MPI=1 -DGMX_VERSION=20250000 \
    -x c "$generated/test-gromacs-version-guard.c.in" >"$tmp/old-gromacs-version.log" 2>&1; then
  echo "GROMACS tMPI adapter unexpectedly accepted a pre-2026 host version" >&2
  exit 1
fi
grep -q 'requires GROMACS 2026.0 or newer' "$tmp/old-gromacs-version.log"

"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -I"$generated" \
  -fsyntax-only -x c++ "$generated/wrapper-header-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -I"$generated" \
  -fsyntax-only -x c++ "$generated/callback-bridges.h"
"${CC:-cc}" -std=c99 -Wall -Wextra -Werror -I"$generated" \
  -fsyntax-only -x c "$generated/provider-callbacks.h"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -I"$tmpi_include" -I"$generated" \
  -fsyntax-only -x c++ "$generated/thread-mpi-symbol-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -I"$tmpi_include" -I"$generated" \
  -fsyntax-only -x c++ "$generated/thread-mpi-entry-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -DPLUMED_MPIW_BUNDLED_THREAD_MPI=1 \
  -I"$repo_root/src/thread_mpi" -I"$generated" \
  -fsyntax-only -x c++ "$generated/thread-mpi-symbol-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -DPLUMED_MPIW_BUNDLED_THREAD_MPI=1 \
  -I"$repo_root/src/thread_mpi" -I"$generated" \
  -fsyntax-only -x c++ "$generated/thread-mpi-entry-check.cpp.in"

MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 -D__PLUMED_HAS_MPI \
  -fsyntax-only -x c++ "$generated/native-header-check.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -fsyntax-only -x c++ "$generated/native-entry-check.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -fsyntax-only -x c++ "$generated/native-signature-check.cpp.in"
mkdir -p "$tmp/mpi1-test"
cp "$generated/test-mpi1-version-header.h.in" "$tmp/mpi1-test/mpi.h"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$tmp/mpi1-test" -iquote "$generated" \
  -fsyntax-only -x c++ "$generated/test-mpi1-version-body.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$tmp/mpi1-test" -iquote "$generated" \
  -fsyntax-only -x c++ "$generated/native-signature-check.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$tmp/mpi1-test" -iquote "$generated" \
  -fsyntax-only -x c++ "$generated/native-entry-check.cpp.in"
OMPI_CXX="${OMPI_CXX:-c++}" "$openmpi_mpicxx" -std=c++17 -D__PLUMED_HAS_MPI \
  -fsyntax-only -x c++ "$generated/native-header-check.cpp.in"
OMPI_CXX="${OMPI_CXX:-c++}" "$openmpi_mpicxx" -std=c++17 \
  -fsyntax-only -x c++ "$generated/native-entry-check.cpp.in"
OMPI_CXX="${OMPI_CXX:-c++}" "$openmpi_mpicxx" -std=c++17 \
  -fsyntax-only -x c++ "$generated/native-signature-check.cpp.in"

echo "generated interface compiled for provider mode, pinned tMPI, MPICH, and Open MPI"
