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
repo_root=$(CDPATH= cd -- "$script_dir/../../.." && pwd)
generated="$repo_root/regtest/basic/rt-mpiwrap-api"
module="$repo_root/src/mpiwrap"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiw-generated.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM
bundled_flag=
if [ "$(CDPATH= cd -- "$tmpi_include" && pwd)" = "$repo_root/src/thread_mpi" ]; then
  bundled_flag=-DPLUMED_MPIW_BUNDLED_THREAD_MPI=1
fi

python3 "$module/generate.py" --check

"${CC:-cc}" -E -I"$generated" -DGMX_THREAD_MPI=1 -DGMX_VERSION=20260000 \
  -x c "$generated/test-gromacs-version-guard.c.in" >/dev/null
if "${CC:-cc}" -E -I"$generated" -DGMX_THREAD_MPI=1 \
    -x c "$generated/test-gromacs-version-guard.c.in" >"$tmp/no-gromacs-version.log" 2>&1; then
  echo "GROMACS tMPI adapter unexpectedly accepted an unknown host version" >&2
  exit 1
fi
grep -q 'requires GMX_VERSION or __PLUMED_WRAPPER_MPI_GROMACS_VERSION' "$tmp/no-gromacs-version.log"
if "${CC:-cc}" -E -I"$generated" -DGMX_THREAD_MPI=1 -DGMX_VERSION=20250000 \
    -x c "$generated/test-gromacs-version-guard.c.in" >"$tmp/old-gromacs-version.log" 2>&1; then
  echo "GROMACS tMPI adapter unexpectedly accepted a pre-2026 host version" >&2
  exit 1
fi
grep -q 'requires GROMACS 2026.0 or newer' "$tmp/old-gromacs-version.log"

"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -I"$generated" \
  -iquote "$module" \
  -fsyntax-only -x c++ "$generated/wrapper-header-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -DPLUMED_TEST_PACKAGE_FIRST -I"$generated" \
  -iquote "$module" \
  -fsyntax-only -x c++ "$generated/wrapper-header-check.cpp.in"
"${CC:-cc}" -std=c99 -Wall -Wextra -Werror -iquote "$module" \
  -fsyntax-only -x c "$module/mpiwrap.h"
"${CXX:-c++}" -std=c++98 -Wall -Wextra -Werror -DPLUMED_MPIW_ABI_ONLY \
  -iquote "$module" -fsyntax-only -x c++ "$module/mpiwrap.h"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -I"$generated" \
  -I"$repo_root/src" -iquote "$module" -fsyntax-only -x c++ "$module/internal.h"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror $bundled_flag -I"$tmpi_include" -I"$generated" \
  -iquote "$module" -I"$repo_root/src" -fsyntax-only -x c++ "$generated/thread-mpi-symbol-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror $bundled_flag -I"$tmpi_include" -I"$generated" \
  -iquote "$module" -I"$repo_root/src" -fsyntax-only -x c++ "$generated/thread-mpi-entry-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -DPLUMED_MPIW_BUNDLED_THREAD_MPI=1 \
  -I"$repo_root/src/thread_mpi" -I"$repo_root/src" -I"$generated" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/thread-mpi-symbol-check.cpp.in"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -DPLUMED_MPIW_BUNDLED_THREAD_MPI=1 \
  -I"$repo_root/src/thread_mpi" -I"$repo_root/src" -I"$generated" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/thread-mpi-entry-check.cpp.in"

MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 -D__PLUMED_HAS_MPI \
  -I"$repo_root/src" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-header-check.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$repo_root/src" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-entry-check.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$repo_root/src" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-signature-check.cpp.in"
mkdir -p "$tmp/mpi1-test"
cp "$generated/test-mpi1-version-header.h.in" "$tmp/mpi1-test/mpi.h"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$tmp/mpi1-test" -I"$repo_root/src" -iquote "$generated" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/test-mpi1-version-body.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$tmp/mpi1-test" -I"$repo_root/src" -iquote "$generated" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-signature-check.cpp.in"
MPICH_CXX="${MPICH_CXX:-c++}" "$mpich_mpicxx" -std=c++17 \
  -I"$tmp/mpi1-test" -I"$repo_root/src" -iquote "$generated" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-entry-check.cpp.in"
OMPI_CXX="${OMPI_CXX:-c++}" "$openmpi_mpicxx" -std=c++17 -D__PLUMED_HAS_MPI \
  -I"$repo_root/src" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-header-check.cpp.in"
OMPI_CXX="${OMPI_CXX:-c++}" "$openmpi_mpicxx" -std=c++17 \
  -I"$repo_root/src" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-entry-check.cpp.in"
OMPI_CXX="${OMPI_CXX:-c++}" "$openmpi_mpicxx" -std=c++17 \
  -I"$repo_root/src" -iquote "$module" \
  -fsyntax-only -x c++ "$generated/native-signature-check.cpp.in"

echo "generated interface compiled for provider mode, pinned tMPI, MPICH, and Open MPI"
