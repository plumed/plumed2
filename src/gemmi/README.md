# Vendored gemmi CIF parser

PLUMED reads mmCIF files (`MOLINFO STRUCTURE=file.cif`, and every other
place where a structure is read with `PDB::read`) with the CIF parser of
[gemmi](https://github.com/project-gemmi/gemmi), release **0.7.5**.
The parser is header-only: `include/gemmi/cif.hpp` and the headers it
includes, among them [PEGTL](https://github.com/taocpp/PEGTL) 2.8.3 as
shipped by gemmi in `include/gemmi/third_party/tao`.

- `include/` is a pristine copy of the release: **never edit these files**.
  Problems are fixed upstream and a new release is imported.
- `./import.sh /path/to/gemmi-X.Y.Z` copies the closure of `cif.hpp`
  (computed by the compiler), `version.hpp` and the licences from an
  extracted release (GitHub tarball or PyPI sdist). The pinned version is
  at the top of the script and is checked against `version.hpp`.
- `cif.h` is PLUMED's only entry point. It includes the copy in `include/`,
  or `<gemmi/cif.hpp>` when PLUMED is configured with `--enable-system-gemmi`
  (a gemmi installed in the system; configure checks that it can be used).
  The upstream headers are treated as system headers, so that PLUMED's
  warning flags (`-Wconversion`, `-Wold-style-cast`, ...) do not apply to them.
- PLUMED does not use gemmi's `Document`: `src/tools/PDBCif.cpp` streams the
  `_atom_site` rows through the grammar `gemmi::cif::rules` with its own
  PEGTL actions, so that a large file is never stored in memory.

Licences: gemmi is MPL-2.0 (`include/gemmi/LICENSE.txt`), PEGTL is MIT
(`include/gemmi/third_party/tao/LICENSE`). Both are compatible with the
LGPL of PLUMED. The files are distributed unmodified.
