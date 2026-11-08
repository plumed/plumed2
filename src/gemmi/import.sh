#! /usr/bin/env bash
# Copy the CIF parser of gemmi (cif.hpp and the headers it includes,
# among them PEGTL) from a gemmi release into include/ of this directory.
# The only change to the files is mechanical and made by this script: the
# namespaces gemmi and tao (PEGTL) are put inside the namespace PLMD, as for
# the other libraries shipped with PLUMED, so that they cannot clash with
# another copy of gemmi or PEGTL linked in the same program (e.g. the one
# bundled in chemfiles, or a host code using gemmi). Never edit the files by
# hand: import a new release instead.
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

# PLMD namespace. Upstream opens and closes the namespaces always in the same way
# (checked below): "namespace gemmi {" ... "} // namespace gemmi" and
# "namespace tao" ... "}  // namespace tao"
for f in $(cd "$tmp/include" && find . -name "*.hpp" | sort); do
  g="$tmp/include/$f"
  open=$(grep -c -e '^namespace gemmi {$' -e '^namespace tao$' "$g" || true)
  close=$(grep -c -e '^} // namespace gemmi$' -e '^}  // namespace tao$' "$g" || true)
  other=$(grep -E '^[[:space:]]*namespace[[:space:]]+(gemmi|tao)([[:space:]]|$|\{)' "$g" | grep -c -v -e '^namespace gemmi {$' -e '^namespace tao$' || true)
  if [[ $open != "$close" || $other != 0 ]]; then
    echo "$f: unexpected namespace layout, the PLMD namespace cannot be added" >&2
    exit 1
  fi
  sed -e 's|^namespace gemmi {$|namespace PLMD { namespace gemmi {|' \
      -e 's|^} // namespace gemmi$|} } // namespace PLMD::gemmi|' \
      -e 's|^namespace tao$|namespace PLMD { namespace tao|' \
      -e 's|^}  // namespace tao$|} }  // namespace PLMD::tao|' "$g" > "$g.tmp"
  mv "$g.tmp" "$g"
done
if grep -rqE -e '(^|[^:A-Za-z_])::(gemmi|tao)::' "$tmp/include"; then
  echo "absolute names ::gemmi:: or ::tao:: found: they would bypass the PLMD namespace" >&2
  exit 1
fi

rm -rf "$here/include"
mv "$tmp/include" "$here/include"
echo "imported gemmi $version: $(find "$here/include" -name "*.hpp" | wc -l) headers into $here/include"
