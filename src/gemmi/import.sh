#! /usr/bin/env bash
# Copy the CIF parser of gemmi (cif.hpp and the headers it includes,
# among them PEGTL) from a gemmi release into include/ of this directory.
# The files are copied unmodified: never edit them, import a new release instead.
#
# usage: ./import.sh /path/to/gemmi-X.Y.Z
# The path is an extracted release (GitHub tarball or the PyPI sdist).
# To change the pinned release, edit "version" below and re-run.
set -euo pipefail

version=0.7.5

if [[ $# != 1 ]]; then
  echo "usage: $0 /path/to/gemmi-$version" >&2
  exit 1
fi

src=$(cd -- "$1" && pwd)
here=$(cd -- "$(dirname -- "$0")" && pwd)

if ! grep -q "#define GEMMI_VERSION \"$version\"" "$src/include/gemmi/version.hpp"; then
  echo "$src is not gemmi $version (see $src/include/gemmi/version.hpp)" >&2
  exit 1
fi

# the closure of cif.hpp, as seen by the compiler
files=$(
  cd "$src/include"
  echo '#include "gemmi/cif.hpp"' |
    ${CXX:-c++} -std=c++17 -MM -I. -x c++ - |
    tr ' \\' '\n\n' | grep '^gemmi/' | sort -u
)
files="$files gemmi/version.hpp"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
for f in $files; do
  mkdir -p "$tmp/include/${f%/*}"
  cp "$src/include/$f" "$tmp/include/$f"
done
cp "$src/LICENSE.txt" "$tmp/include/gemmi/LICENSE.txt"
cp "$src/include/gemmi/third_party/tao/LICENSE" "$tmp/include/gemmi/third_party/tao/LICENSE"
cp "$src/include/gemmi/third_party/tao/NOTES" "$tmp/include/gemmi/third_party/tao/NOTES"

rm -rf "$here/include"
mv "$tmp/include" "$here/include"
echo "imported gemmi $version: $(find "$here/include" -name "*.hpp" | wc -l) headers into $here/include"
