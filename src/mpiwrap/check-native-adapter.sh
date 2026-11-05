#!/bin/sh
set -eu

module_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
mpich_cc=${MPICH_CC_WRAPPER:?set MPICH_CC_WRAPPER to an absolute MPI compiler path}
mpich_cxx=${MPICH_CXX_WRAPPER:?set MPICH_CXX_WRAPPER to an absolute MPI compiler path}
ompi_cc=${OMPI_CC_WRAPPER:?set OMPI_CC_WRAPPER to an absolute MPI compiler path}
ompi_cxx=${OMPI_CXX_WRAPPER:?set OMPI_CXX_WRAPPER to an absolute MPI compiler path}
clang_c=${CLANG_C:-/usr/bin/clang}
clang_cxx=${CLANG_CXX:-/usr/bin/clang++}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiw-native-adapter.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

for provider in mpich openmpi; do
  case $provider in
    mpich) cc=$mpich_cc; cxx=$mpich_cxx; export MPICH_CC=$clang_c MPICH_CXX=$clang_cxx ;;
    openmpi) cc=$ompi_cc; cxx=$ompi_cxx; export OMPI_CC=$clang_c OMPI_CXX=$clang_cxx ;;
  esac
  extra_cflags=
  if [ "$provider" = mpich ]; then extra_cflags=-DPLUMED_MPIW_EXPECT_NONNULL_STATUS_IGNORE=1; fi
  for capacity in default 16; do
    capacity_flag=
    if [ "$capacity" = 16 ]; then capacity_flag=-DPLUMED_MPIW_HANDLE_CAPACITY=16; fi
    "$cc" -std=c99 -Wall -Wextra -Werror $extra_cflags $capacity_flag -x c \
      "$module_dir/test-native-adapter.c.in" -o "$tmp/test-$provider-$capacity-c"
    "$tmp/test-$provider-$capacity-c"
    "$cc" -x c -std=c99 -Wall -Wextra -Werror $capacity_flag -DPLUMED_MPIW_NATIVE_ADAPTER=1 \
      -c \
      "$module_dir/mpiwrap_implementation.cpp" -o "$tmp/adapter-$provider-$capacity-c99.o"
    "$cxx" -std=c++98 -Wall -Wextra -Werror $capacity_flag -DPLUMED_MPIW_NATIVE_ADAPTER=1 \
      -c \
      "$module_dir/mpiwrap_implementation.cpp" -o "$tmp/adapter-$provider-$capacity-cxx98.o"
    "$cc" -x c -std=c99 -Wall -Wextra -Werror -DPLUMED_MPIW_HANDLE_CAPACITY=1 \
      -DPLUMED_MPIW_STATUS_CAPACITY=1 -DPLUMED_MPIW_NATIVE_ADAPTER=1 \
      -c \
      "$module_dir/mpiwrap_implementation.cpp" -o "$tmp/adapter-$provider-tiny-wrapper-c99.o"
    "$cxx" -std=c++98 -Wall -Wextra -Werror -DPLUMED_MPIW_HANDLE_CAPACITY=1 \
      -DPLUMED_MPIW_STATUS_CAPACITY=1 -DPLUMED_MPIW_NATIVE_ADAPTER=1 \
      -c \
      "$module_dir/mpiwrap_implementation.cpp" -o "$tmp/adapter-$provider-tiny-wrapper-cxx98.o"
    for host_language in c cxx98; do
      case $host_language in
        c) host_compiler=$cc; host_standard=c99; host_input_language=c ;;
        cxx98) host_compiler=$cxx; host_standard=c++98; host_input_language=c++ ;;
      esac
      "$host_compiler" -std="$host_standard" -Wall -Wextra -Werror $capacity_flag \
        -D__PLUMED_WRAPPER_IMPLEMENTATION=1 -DPLUMED_MPIW_HOST_ADAPTER=1 \
        -D__PLUMED_HAS_DLOPEN=1 -iquote "$module_dir" -x "$host_input_language" -c \
        "$module_dir/test-host-header-adapter.c.in" \
        -o "$tmp/host-$provider-$capacity-$host_language.o"
    done
  done
  for language in c99 cxx98; do
    # Both bodies implement the C API only, including the C++98 object.
    # Linking with the C driver must not require legacy MPI C++ bindings.
    linker=$cc
    library="$tmp/provider-$provider-$language.dylib"
    "$linker" -shared "$tmp/adapter-$provider-default-$language.o" -o "$library"
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

if [ -n "${GMX_THREAD_MPI_INCLUDE:-}" ]; then
  for host_language in c cxx98; do
    case $host_language in
      c) host_compiler=$clang_c; host_standard=c99; host_input_language=c ;;
      cxx98) host_compiler=$clang_cxx; host_standard=c++98; host_input_language=c++ ;;
    esac
    "$host_compiler" -std="$host_standard" -Wall -Wextra -Werror \
      -DGMX_THREAD_MPI=1 -DPLUMED_MPIW_GROMACS_VERSION=20260000 \
      -DPLUMED_MPIW_TEST_THREAD_MPI=1 \
      -DSIZEOF_LONG_LONG_INT=8 -D__PLUMED_WRAPPER_IMPLEMENTATION=1 \
      -DPLUMED_MPIW_HOST_ADAPTER=1 -D__PLUMED_HAS_DLOPEN=1 \
      -I"$GMX_THREAD_MPI_INCLUDE" -iquote "$module_dir" -x "$host_input_language" -c \
      "$module_dir/test-host-header-adapter.c.in" \
      -o "$tmp/host-thread-mpi-$host_language.o"
  done
fi
