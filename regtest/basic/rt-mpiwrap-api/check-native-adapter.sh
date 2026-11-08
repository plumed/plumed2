#!/bin/sh
set -eu

fixture_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$fixture_dir/../../.." && pwd)
module_dir="$root/src/mpiwrap"
mpich_cc=${MPICH_CC_WRAPPER:?set MPICH_CC_WRAPPER to an absolute MPI compiler path}
mpich_cxx=${MPICH_CXX_WRAPPER:?set MPICH_CXX_WRAPPER to an absolute MPI compiler path}
ompi_cc=${OMPI_CC_WRAPPER:?set OMPI_CC_WRAPPER to an absolute MPI compiler path}
ompi_cxx=${OMPI_CXX_WRAPPER:?set OMPI_CXX_WRAPPER to an absolute MPI compiler path}
clang_c=${CLANG_C:-/usr/bin/clang}
clang_cxx=${CLANG_CXX:-/usr/bin/clang++}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiw-native-adapter.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM
case $(uname -s) in
  Darwin) shared_flag=-dynamiclib ;;
  *) shared_flag=-shared ;;
esac

for provider in mpich openmpi; do
  case $provider in
    mpich) cc=$mpich_cc; cxx=$mpich_cxx; export MPICH_CC=$clang_c MPICH_CXX=$clang_cxx ;;
    openmpi) cc=$ompi_cc; cxx=$ompi_cxx; export OMPI_CC=$clang_c OMPI_CXX=$clang_cxx ;;
  esac
  MPICC="$cc" MPICXX="$cxx" bash "$fixture_dir/check-host-header-linkage.sh"
  extra_cflags=
  if [ "$provider" = mpich ]; then extra_cflags=-DPLUMED_MPIW_EXPECT_NONNULL_STATUS_IGNORE=1; fi
  for capacity in default 16; do
    capacity_flag=
    if [ "$capacity" = 16 ]; then capacity_flag=-DPLUMED_MPIW_HANDLE_CAPACITY=16; fi
    "$cc" -std=c99 -Wall -Wextra -Werror $extra_cflags $capacity_flag -iquote "$module_dir" -x c \
      "$fixture_dir/test-native-adapter.c.in" -o "$tmp/test-$provider-$capacity-c"
    "$tmp/test-$provider-$capacity-c"
    "$cxx" -std=c++98 -Wall -Wextra -Werror $capacity_flag -iquote "$module_dir" -x c++ \
      "$fixture_dir/test-native-adapter.c.in" -o "$tmp/test-$provider-$capacity-cxx98"
    "$tmp/test-$provider-$capacity-cxx98"
    for host_language in c cxx98; do
      case $host_language in
        c) host_compiler=$cc; host_standard=c99; host_input_language=c ;;
        cxx98) host_compiler=$cxx; host_standard=c++98; host_input_language=c++ ;;
      esac
      "$host_compiler" -std="$host_standard" -Wall -Wextra -Werror $capacity_flag \
        -D__PLUMED_WRAPPER_IMPLEMENTATION=1 -D__PLUMED_WRAPPER_MPI_HOST_ADAPTER=1 \
        -D__PLUMED_HAS_DLOPEN=1 -iquote "$module_dir" -x "$host_input_language" -c \
        "$fixture_dir/test-host-header-adapter.c.in" \
        -o "$tmp/host-$provider-$capacity-$host_language.o"
    done
  done
  for anonymous in 0 1; do
    for language in c99 cxx98; do
      # The public wrapper embeds a C-only adapter body in both language modes.
      linker=$cc
      library="$tmp/provider-$provider-$language.dylib"
      case $language in
        c99) std=c99; input_language=c ;;
        cxx98) std=c++98; input_language=c++ ;;
      esac
      "$linker" "$shared_flag" -fPIC -std="$std" -iquote "$module_dir" \
        -D__PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE=$anonymous \
        -x "$input_language" "$fixture_dir/test-native-adapter.c.in" -o "$library"
      python3 - "$library" <<'PYTHON'
import ctypes
import os
import sys
ctypes.CDLL(sys.argv[1], mode=os.RTLD_NOW)
PYTHON
      if ! nm -g "$library" | grep -E '[[:space:]]_?plumed_mpiw_provider_factory$' >/dev/null; then
        echo "$provider $language provider is missing its unmangled factory" >&2
        exit 1
      fi
      if nm -u "$library" | grep -E '[[:space:]]_?plumed([[:alnum:]_]*)([[:space:]]|$)' >/dev/null; then
        echo "$provider $language provider has an unresolved PLUMED dependency" >&2
        exit 1
      fi
      if command -v otool >/dev/null 2>&1 && otool -L "$library" | grep -i plumed >/dev/null; then
        echo "$provider $language provider links a PLUMED library" >&2
        exit 1
      fi
    done
  done
done

if [ -n "${GMX_THREAD_MPI_INCLUDE:-}" ]; then
  for long_long_support in enabled disabled; do
    for host_language in c cxx98; do
      case $host_language in
        c) host_compiler=$clang_c; host_standard=c99; host_input_language=c ;;
        cxx98) host_compiler=$clang_cxx; host_standard=c++98; host_input_language=c++ ;;
      esac
      long_long_flag=
      if [ "$long_long_support" = enabled ]; then
        long_long_flag=-DSIZEOF_LONG_LONG_INT=8
      fi
      "$host_compiler" -std="$host_standard" -Wall -Wextra -Werror \
        -DGMX_THREAD_MPI=1 -D__PLUMED_WRAPPER_MPI_GROMACS_VERSION=20260000 \
        -DPLUMED_MPIW_TEST_THREAD_MPI=1 \
        $long_long_flag -D__PLUMED_WRAPPER_IMPLEMENTATION=1 \
        -D__PLUMED_WRAPPER_MPI_HOST_ADAPTER=1 -D__PLUMED_HAS_DLOPEN=1 \
        -I"$GMX_THREAD_MPI_INCLUDE" -iquote "$module_dir" -x "$host_input_language" -c \
        "$fixture_dir/test-host-header-adapter.c.in" \
        -o "$tmp/host-thread-mpi-$long_long_support-$host_language.o"
    done
  done
fi
