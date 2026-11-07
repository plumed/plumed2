#!/usr/bin/env bash
set -euo pipefail

# Use the installed public interface. No PLUMED source rebuild or MPI-dependent
# kernel is needed; the two variants differ only in these small binaries.
export MPICH_CC="$CC" OMPI_CC="$CC"
flags=()
read -r -a compile_flags <<< "${CPPFLAGS:-} ${CFLAGS:-}"
read -r -a link_flags <<< "${LDFLAGS:-}"
for flag in ${compile_flags[@]+"${compile_flags[@]}"}; do flags+=(--cflag "$flag"); done
for flag in ${link_flags[@]+"${link_flags[@]}"}; do flags+=(--ldflag "$flag"); done
plumed mkmpi ${flags[@]+"${flags[@]}"} --compiler "$PREFIX/bin/mpicc" \
  --out "$PREFIX/lib/libplumed-mpi-provider$SHLIB_EXT"
plumed mkmpi ${flags[@]+"${flags[@]}"} --executable --link shared --compiler "$PREFIX/bin/mpicc" \
  --out "$PREFIX/bin/plumed-mpi"
