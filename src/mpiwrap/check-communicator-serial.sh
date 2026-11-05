#!/usr/bin/env bash
set -euo pipefail

if [[ $# == 1 ]]; then
  source_root=$(cd "$(dirname "$0")/../.." && pwd)
  kernel_path=$1
elif [[ $# == 2 ]]; then
  source_root=$(cd "$1" && pwd)
  kernel_path=$2
else
  echo "usage: $0 [/path/to/source-root] /path/to/libplumedKernel.dylib-or-so" >&2
  exit 2
fi

module_dir=$(cd "$(dirname "$0")" && pwd)
if [[ ! -f $kernel_path || ! -f $source_root/src/tools/Communicator.h ]]; then
  echo "source root and kernel shared library must exist" >&2
  exit 2
fi
kernel=$(cd "$(dirname "$kernel_path")" && pwd)/$(basename "$kernel_path")
kernel_dir=$(dirname "$kernel")
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/mpiwrap-communicator-serial.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM

"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -Wno-unused-parameter \
  -D__PLUMED_MPIWRAP=1 -I"$source_root/src" \
  -x c++ "$module_dir/test-communicator-serial.cpp.in" -x none "$kernel" \
  "-Wl,-rpath,$kernel_dir" -o "$tmpdir/test-communicator-serial"
"$tmpdir/test-communicator-serial"

echo "serial Communicator fallback stayed nonselecting and preserved local behavior"
