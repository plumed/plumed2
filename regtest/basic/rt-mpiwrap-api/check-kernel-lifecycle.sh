#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
CXX=${CXX:-c++}
TMP=$(mktemp -d "${TMPDIR:-/tmp}/mpiwrap-kernel-test.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

case "$(uname -s)" in
  Darwin) SHARED_FLAGS=(-dynamiclib) ;;
  *) SHARED_FLAGS=(-shared) ;;
esac

"$CXX" -std=c++17 -fPIC "${SHARED_FLAGS[@]}" \
  -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/fake-kernel.cpp.in" -o "$TMP/libkernel-api13.so"
"$CXX" -std=c++17 -fPIC "${SHARED_FLAGS[@]}" \
  -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/fake-old-kernel.cpp.in" -o "$TMP/libkernel-api10.so"
"$CXX" -std=c++17 -fPIC -DTEST_REJECT_INSTALL "${SHARED_FLAGS[@]}" \
  -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/fake-kernel.cpp.in" -o "$TMP/libkernel-reject.so"
"$CXX" -std=c++17 -D__PLUMED_HAS_DLOPEN -I"$ROOT/src" \
  -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/test-kernel-lifecycle.cpp.in" -ldl -o "$TMP/test-kernel-lifecycle"
"$TMP/test-kernel-lifecycle" "$TMP/libkernel-api13.so" \
  "$TMP/libkernel-api10.so" "$TMP/libkernel-reject.so" \
  "$TMP/kernel-unloaded" "$TMP/object-created"
"$CXX" -std=c++17 -D__PLUMED_HAS_DLOPEN -I"$ROOT/src" \
  -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/test-forward-provider.cpp.in" -ldl -o "$TMP/test-forward-provider"
"$TMP/test-forward-provider" "$TMP/libkernel-api13.so" "$TMP/libkernel-api10.so" \
  "$TMP/libkernel-reject.so"
# Upstream APIs 11 and 12 add unrelated commands, not MPI injection.
for api in 11 12; do
  "$CXX" -std=c++17 -fPIC -DTEST_API_VERSION="$api" "${SHARED_FLAGS[@]}" \
    -x c++ "$ROOT/regtest/basic/rt-mpiwrap-api/fake-old-kernel.cpp.in" -o "$TMP/libkernel-api$api.so"
  "$TMP/test-forward-provider" "$TMP/libkernel-api13.so" "$TMP/libkernel-api$api.so" \
    "$TMP/libkernel-reject.so"
done
echo "kernel create/inject/release lifecycle and API 10/11/12 compatibility passed"
