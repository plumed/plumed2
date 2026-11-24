#!/usr/bin/env bash
set -euo pipefail
if (($# != 1)) || [[ ! -f $1 ]]; then
  echo "usage: $0 /path/to/configured/libplumed.a" >&2
  exit 2
fi
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-thread-mpi-lifecycle.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
CXX=${CXX:-clang++}
cp "$ROOT/regtest/basic/rt-mpiwrap-launcher/test-thread-mpi-lifecycle.cpp.in" "$TASK_DIR/thread-mpi-lifecycle.cpp"
"$CXX" -std=c++17 -D__PLUMED_MPIWRAP=1 -D__PLUMED_HAS_THREAD_MPI=1 \
  -I"$ROOT/src" -I"$ROOT/src/config" \
  "$TASK_DIR/thread-mpi-lifecycle.cpp" "$1" \
  -lz -ldl -llapack -lblas -pthread -o "$TASK_DIR/thread-mpi-lifecycle"
cat > "$TASK_DIR/plumed.dat" <<'PLUMED'
UNITS LENGTH=A ENERGY=kcal/mol TIME=ps
PLUMED
case "$(uname -s)" in
  Darwin) compatibility_shared=-dynamiclib ;;
  *) compatibility_shared=-shared ;;
esac
for api in 10 11 12; do
  "$CXX" -std=c++17 -fPIC "$compatibility_shared" -DTEST_API_VERSION="$api" \
    -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/fake-old-kernel.cpp.in" \
    -o "$TASK_DIR/old-api$api.so"
done
(cd "$TASK_DIR" && "$TASK_DIR/thread-mpi-lifecycle" \
  "$TASK_DIR/old-api10.so" "$TASK_DIR/old-api11.so" "$TASK_DIR/old-api12.so")
echo "thread-MPI sentinel guard and repeated-launch lifecycle checks passed"
