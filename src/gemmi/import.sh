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
# The path is an extracted release (GitHub tarball or the PyPI sdist). Needs python3.
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

# the closure of cif.hpp: every header included with #include "...", also inside #if
# blocks, so that the copy works on every platform (e.g. utf.hpp and demangle_nop.hpp
# are used only on Windows/MSVC). The PEGTL headers that gemmi does not ship
# (see third_party/tao/NOTES, e.g. file_mapper_win32.hpp: no mmap on Windows) are skipped.
files=$(
  cd "$src/include"
  python3 - <<'PYEOF'
import os, re, sys
todo, seen = ["gemmi/cif.hpp"], set()
inc = re.compile(r'^\s*#\s*include\s*"([^"]+)"')
while todo:
    f = todo.pop()
    if f in seen:
        continue
    seen.add(f)
    for line in open(f, encoding="utf-8"):
        m = inc.match(line)
        if m:
            g = os.path.normpath(os.path.join(os.path.dirname(f), m.group(1)))
            if os.path.exists(g):
                todo.append(g)
            else:
                print("not in the release, skipped: " + g + " (included by " + f + ")", file=sys.stderr)
print("\n".join(sorted(seen)))
PYEOF
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
