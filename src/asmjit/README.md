# Vendored AsmJit 1.9.0

Imported from `libraries/asmjit` in OpenMM 8.6.1 (released 2026-09-11), commit
`b399af4725573963b46d6c1083fdcf7a37615857`.
This is the version paired with OpenMM's Lepton ARM64 and x86 JIT backends.

Check out OpenMM tag `8.6.1` before importing.

```sh
./import.sh /path/to/openmm/libraries/asmjit
# A standalone checkout of the matching AsmJit 1.9.0 sources also works:
./import.sh /path/to/asmjit
```

The importer requires Python 3. It checks the upstream version and duplicate
filenames before replacing sources, flattens directories and quoted includes,
renames `string.h` to `moved_string.h` to avoid shadowing the C library header,
and isolates all AsmJit namespaces and qualified references under `PLMD`.
Sources compile only when `__PLUMED_HAS_ASMJIT` is enabled. Static embedding is
selected, upstream licenses are retained, and Clang warning suppressions are
scoped using upstream's balanced diagnostic push/pop macros.

For a reproduction check without modifying this directory:

```sh
python3 import.py /path/to/openmm/libraries/asmjit --output-dir /tmp/asmjit-reimport
```

Compare the generated `.h` and `.cpp` files byte for byte. The import does not
run `header.sh` or stage changes in Git. A newer AsmJit API requires a separate
port of both the importer and Lepton's JIT implementation.
