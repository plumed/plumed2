# Embedded thread-MPI

This module contains the namespaced GROMACS thread-MPI import. [`UPSTREAM.json`](UPSTREAM.json)
is the single provenance record: it pins the resolved source commit and records
SHA256 checksums for each original upstream file before import transformations.
The license is copied from the pinned upstream tree. Like other vendored
libraries, this module is excluded from the PLUMED-owned compiler-warning
gate and coverage metrics. The `mpiwrap` adapter remains subject to both. The imported runtime has no
PLUMED module dependencies (`USE=`); `mpiwrap` depends on it, rather than
the runtime depending on `core`.

`reimport.py` discovers implementation sources and headers from the pinned
upstream directory, flattens them into this directory, rewrites internal
includes, and places C++ declarations and definitions in `PLMD::thread_mpi`.
When flattened basenames collide, it derives stable names from their source
roles and directory paths. Imported translation units retain the
`__PLUMED_HAS_THREAD_MPI` guard. The importer formats generated C++ and headers
with the project's required `astyle` options.

To verify the recorded import against a local GROMACS checkout:

```sh
python3 src/thread_mpi/reimport.py /path/to/gromacs --check
```

To update from a local GROMACS tag or commit, pass it directly. The script
resolves the commit, discovers source files, refreshes checksums, and writes
the provenance record; no manual edits to `UPSTREAM.json` are needed. This
also bootstraps the manifest if it is missing.

```sh
python3 src/thread_mpi/reimport.py /path/to/gromacs --ref v2026.4
```

To check a candidate ref without changing the import or its provenance, add
`--check`:

```sh
python3 src/thread_mpi/reimport.py /path/to/gromacs --ref v2026.4 --check
```

Review any build or source-layout changes in the upstream release.

The standalone runtime tests live in
[`regtest/basic/rt-mpiwrap-launcher`](../../regtest/basic/rt-mpiwrap-launcher/README.md).
Build them from the repository root with:

```sh
cmake -S regtest/basic/rt-mpiwrap-launcher -B /tmp/thread-mpi-build
cmake --build /tmp/thread-mpi-build
ctest --test-dir /tmp/thread-mpi-build --output-on-failure
python3 regtest/basic/rt-mpiwrap-launcher/check-symbol-isolation.py /path/to/gromacs
```

MPI provider design, wrapper behavior, and validation limits are documented in
[`../mpiwrap/README.md`](../mpiwrap/README.md).

The standalone suite includes a GROMACS-style `MPI_*` macro-surface test using
`mpi_bindings.h`: 1, 2 and 4 ranks exchange a committed contiguous datatype and
check communicator and predefined-value bindings. This uses the embedded tMPI
sources; it does not build GROMACS or claim to test GROMACS integration.
Sequential internal launches with different rank counts are also checked.
