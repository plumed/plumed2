#!/usr/bin/env bash
# Validate injection priority and the explicit troubleshooting override.
set -euo pipefail
if (($# != 2)); then
  echo 'usage: check-mpi-provider-force.sh KERNEL PROVIDER' >&2
  exit 2
fi
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-provider-force.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
case "$(uname -s)" in
  Darwin) shared_flags=(-dynamiclib); loader_flags=() ;;
  *) shared_flags=(-shared); loader_flags=(-ldl) ;;
esac
"${CXX:-c++}" -std=c++17 -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/test-provider-force.cpp.in" \
  ${loader_flags[@]+"${loader_flags[@]}"} -o "$TASK_DIR/probe"
cat > "$TASK_DIR/bypassed.cpp" <<'SOURCE'
#include <cstdio>
#include <cstdlib>
struct plumed_mpi_provider;
namespace PLMD {
extern "C" const plumed_mpi_provider *plumed_mpiw_provider_factory() { return 0; }
__attribute__((destructor)) void bypassedProviderUnloaded() {
  const char *path=std::getenv("MPIW_BYPASSED_MARKER");
  if(path) {
    std::FILE *file=std::fopen(path,"w");
    if(file) std::fclose(file);
  }
}
}
SOURCE
"${CXX:-c++}" -fPIC "${shared_flags[@]}" "$TASK_DIR/bypassed.cpp" -o "$TASK_DIR/bypassed.so"
# Normal injection wins over an unusable environment selection.
PLUMED_MPI_WRAPPER="$TASK_DIR/missing.so" PLUMED_MPI_WRAPPER_FORCE=0 \
  "$TASK_DIR/probe" "$1" "$2"
# Force bypasses even an invalid incoming descriptor. The host releases its
# bypassed-library reference after destroying the PLUMED object and kernel.
MPIW_BYPASSED_MARKER="$TASK_DIR/unloaded" PLUMED_MPI_WRAPPER="$2" \
  PLUMED_MPI_WRAPPER_FORCE=1 "$TASK_DIR/probe" "$1" "$TASK_DIR/bypassed.so"
test -f "$TASK_DIR/unloaded"
for setting in missing-path bad-force invalid-provider; do
  case "$setting" in
    missing-path) selected=; force=1 ;;
    bad-force) selected=$2; force=invalid ;;
    invalid-provider) selected=$TASK_DIR/missing.so; force=1 ;;
  esac
  if PLUMED_MPI_WRAPPER="$selected" PLUMED_MPI_WRAPPER_FORCE="$force" \
      "$TASK_DIR/probe" "$1" "$2" > "$TASK_DIR/error.log" 2>&1; then
    echo "invalid forced provider succeeded: $setting" >&2
    exit 1
  fi
  grep -E 'PLUMED_MPI_WRAPPER|MPI provider' "$TASK_DIR/error.log" >/dev/null
done
echo 'injection priority, forced selection, and bypassed-library ownership passed'
