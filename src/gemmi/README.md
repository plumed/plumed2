# Vendored gemmi CIF parser

PLUMED reads mmCIF files (`MOLINFO STRUCTURE=file.cif`, and every other
place where a structure is read with `PDB::read`) with the CIF parser of
[gemmi](https://github.com/project-gemmi/gemmi), release **0.7.5**.
The parser is header-only: `include/gemmi/cif.hpp` and the headers it
includes, among them [PEGTL](https://github.com/taocpp/PEGTL) 2.8.3 as
shipped by gemmi in `include/gemmi/third_party/tao`.

- `include/` is a copy of the release with one mechanical change made by
  `import.sh`: the namespaces `gemmi` and `tao` (PEGTL) are put inside the
  namespace `PLMD`, as for the other libraries shipped with PLUMED. The copy
  thus cannot clash with another gemmi or PEGTL linked in the same program
  (chemfiles bundles its own gemmi; a host code or a python module may use
  gemmi too). **Never edit these files by hand**: problems are fixed upstream
  and a new release is imported.
- `./import.sh /path/to/gemmi-X.Y.Z` copies the closure of `cif.hpp`
  (every header reached through `#include "..."`, also inside `#if` blocks, so
  that the copy works on every platform: e.g. `utf.hpp` and `demangle_nop.hpp`
  are used only on Windows/MSVC; the few PEGTL headers that gemmi does not ship
  are skipped), `version.hpp` and the licences from an
  extracted release (GitHub tarball or PyPI sdist). The pinned version is
  at the top of the script and is checked against `version.hpp`.
- `cif.h` is PLUMED's only entry point. It includes the copy in `include/`,
  or `<gemmi/cif.hpp>` when PLUMED is configured with `--enable-system-gemmi`
  (a gemmi installed in the system; configure checks that it can be used).
  A system gemmi stays in the global namespaces, and `cif.h` makes
  `PLMD::gemmi` and `PLMD::tao` aliases of them, so the code is the same.
  The upstream headers are treated as system headers, so that PLUMED's
  warning flags (`-Wconversion`, `-Wold-style-cast`, ...) do not apply to them.
- PLUMED does not use gemmi's `Document`: `src/tools/PDBCif.cpp` streams the
  `_atom_site` rows through the grammar `gemmi::cif::rules` with its own
  PEGTL actions, so that a large file is never stored in memory.

Licences: gemmi is MPL-2.0 (`include/gemmi/LICENSE.txt`), PEGTL is MIT
(`include/gemmi/third_party/tao/LICENSE`). Both are compatible with the
LGPL of PLUMED. Apart from the namespace lines changed by `import.sh`,
the files are distributed unmodified.
