#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_root=$(CDPATH= cd -- "$test_dir/../../.." && pwd)
build_dir="$source_root/.tmp/wrapper-shapes-check"
compiler=${CXX:-c++}
mkdir -p "$build_dir/include"
# Keep the check independent of configure-generated include symlinks and of
# any installed kernel's source directory.
ln -sfn "$source_root/src" "$build_dir/include/plumed"

for language in c++11 c++14 c++17; do
  for strict in 0 1; do
    executable="$build_dir/check-$language-$strict"
    "$compiler" -std="$language" -pedantic-errors -Wall -Wextra -Werror \
      -Wno-unknown-pragmas -Wno-unused-parameter -D__PLUMED_WRAPPER_CXX_DETECT_SHAPES_STRICT="$strict" \
      -I "$build_dir/include" "$test_dir/main.cpp" -o "$executable"
    (cd "$build_dir" && "$executable")
  done
done

for rejection in nested bool nonpointer volatile; do
  source="$test_dir/reject-$rejection.cpp.in"
  log="$build_dir/reject-$rejection.log"
  if "$compiler" -std=c++11 -pedantic-errors -I "$build_dir/include" \
      -x c++ -c "$source" -o "$build_dir/reject-$rejection.o" >"$log" 2>&1; then
    echo "unexpectedly accepted $source" >&2
    exit 1
  fi
  case "$rejection" in
    nested) diagnostic="cmd cannot flatten nested dynamic containers" ;;
    bool|volatile) diagnostic="no matching.*SafePtr" ;;
    nonpointer) diagnostic="no matching.*cmd_helper_with_shape" ;;
  esac
  if ! grep -E "$diagnostic" "$log" >/dev/null; then
    cat "$log" >&2
    exit 1
  fi
done
