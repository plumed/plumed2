# MPI wrapper and embedded thread-MPI tests

This directory contains the MPI wrapper regtest and a standalone CMake suite
for PLUMED's embedded thread-MPI runtime. The imported GROMACS runtime remains
in [`src/thread_mpi`](../../../src/thread_mpi); its provenance and reimport
workflow are documented there.

Configure, build, and run the suite from the repository root:

```sh
cmake -S regtest/basic/rt-mpiwrap-launcher -B /tmp/thread-mpi-build
cmake --build /tmp/thread-mpi-build
ctest --test-dir /tmp/thread-mpi-build --output-on-failure
```

CMake registers 15 tests: launch, MPI bindings, and internal-launcher checks at
1, 2, and 4 ranks; callback-exception and unlaunched-serial checks; runtime
dispatch at 1, 2, and 4 ranks; and an injected-provider smoke test across a
host/kernel DSO boundary. The fixture sources stay here. CMake compiles the
imported runtime from `src/thread_mpi` and the wrapper sources from
`src/mpiwrap`.

To check that imported PLUMED symbols remain isolated from upstream GROMACS
tMPI, run:

```sh
python3 regtest/basic/rt-mpiwrap-launcher/check-symbol-isolation.py /path/to/gromacs
```
