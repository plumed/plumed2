# Conda packaging

These are PLUMED's own CI/channel recipes. They prototype optional extension
packages; they do not change conda-forge's recipes or existing package names.
Runtime selection, host compatibility, lifecycle rules and provider limits are
covered by the [MPI wrapper documentation](../src/mpiwrap/README.md).

| Package | Contents | External MPI dependency |
| --- | --- | --- |
| `plumed` | One serial-built kernel, wrappers, headers, CLI, and embedded thread-MPI | None |
| `plumed-mpi` (`mpi_mpich_*` variant) | MPICH provider library and `plumed-mpi` executable | MPICH |
| `plumed-mpi` (`mpi_openmpi_*` variant) | Open MPI provider library and `plumed-mpi` executable | Open MPI |
| `py-plumed` | Python wrapper, depending on the same `plumed` | None |

`plumed` retains its name; there is no `plumed-kernel` package or `nompi`
add-on. It is built once per platform, explicitly disabling native MPI and
using the non-MPI FFTW variant. Adding/removing an extension must not replace
its kernel. Optional source modules are delivered as `plumed-ext-<module>` packages,
without an MPI variant dimension. See [extension builds and loading](../extensions/README.md).

Both MPI variants depend on the matching PLUMED version, use the installed
`Plumed.h`/standard main via `plumed mkmpi`, and install only:

- `lib/libplumed-mpi-provider.so` (`.dylib` on macOS);
- `bin/plumed-mpi`.

The two providers are variants of one package, so only one is installed at a
time. MPI variants follow the [conda-forge MPI recipe conventions](https://conda-forge.org/docs/maintainer/knowledge_base/#message-passing-interface-mpi).
The base recipe has no MPI variant dimension and the Python recipe does not
multiply its builds by MPI implementation.

## Usage

Select this CI/channel's packages when trying the prototype:

```sh
conda install -c CHANNEL plumed
# Optional: choose a provider to match the rest of the environment.
conda install -c CHANNEL plumed-mpi mpich
# Or select Open MPI instead:
conda install -c CHANNEL plumed-mpi openmpi
mpiexec -n 4 plumed-mpi driver --plumed plumed.dat --ixyz trajectory.xyz
```

The ordinary `plumed` executable keeps its serial/internal-thread launcher.
`plumed-mpi` is the standard PLUMED main rebuilt with MPI: it owns MPI startup
and injects the implementation. Hosts using older wrappers can instead pass
their native communicator; the kernel then discovers the prebuilt provider
beside itself. Neither route needs `PLUMED_MPI_WRAPPER` or a runtime compiler.
The matching MPI runtime still belongs to the host environment.

A provider file is not selected merely by constructing an ordinary PLUMED
object. Internal thread-MPI remains available through its explicit launcher.
Removing `plumed-mpi` removes its executable/provider while leaving the base
installation usable. This prototype does not make the serial `plumed`
executable initialize external MPI automatically.

## Building and checking

In a conda environment with `conda-build`, activate it and run:

```sh
export VERSION=none  # CI development builds; tags use their release version
make -C conda
make -C conda test
```

The build script builds `plumed` once, then both small `plumed-mpi` variants,
then the existing Python variants. Local packages take priority when building
and testing dependents. No packages are uploaded by these commands; the
existing upload target remains a separate, token-controlled release step.

Each MPI package's relocated test environment compiles a small native host
without injection. Both this host and the installed injecting executable run
a two-rank distance calculation. MPI compiler selection is then deliberately
set to `false`, so missing automatic discovery cannot be hidden by compilation
of another adapter. The same input also runs with explicit internal thread-MPI
while the provider file is installed. Compiler packages are test/build
requirements, not runtime requirements of `plumed-mpi`.

`make -C conda test` creates fresh environments from the local channel. For
each implementation it checks that the base has no external MPI dependency,
adds `plumed-mpi` plus the selected runtime, exercises the MPI executable,
removes the add-on, and checks that the kernel bytes stayed unchanged and
ordinary/internal-thread capabilities remain available. CI runs these checks
on Linux x86_64, native Linux aarch64 (`ubuntu-24.04-arm`), and macOS.
The workflow uses `conda-ci`'s `master` branch for the bootstrap on all platforms.

On macOS, the recipes mark the main executables as `no_link`: conda installs
copies instead of hardlinks to its package cache. Otherwise dyld can resolve
relative library paths through another hardlink, loading a kernel from the cache
instead of the selected environment. Library installation is unchanged.
