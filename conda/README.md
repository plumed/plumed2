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
| `plumed-ext-<module>` | Separately loaded optional module and registration | None |
| `py-plumed` | Python wrapper, depending on the same `plumed` | None |
| `plumed-mklib` | Runtime C++ compiler and compilation scripts for the exact kernel | None |

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

The build script builds `plumed` once, generates the extension inventory and
recipes, then builds matching Python/tooling, optional module, MPI provider and
plugin variants, including opted-in empty compatibility packages. Module-local `conda.mk` files declare
native conda requirements; licenses come from module license files. Optional
libraries that are unavailable do not exclude modules during this proof of concept.
To invoke conda-build directly on extensions, first run
`python extensions/conda_recipe.py` from the repository root. Local packages take priority when building
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

## Optional runtime compilation

The main `plumed` package has no runtime compiler dependency. Install
`plumed-mklib` to use `plumed mklib` or `LOAD FILE=action.cpp`. It depends on
that exact kernel build and a compatible C++ compiler. Activate the environment
before compiling; `PLUMED_MKLIB_CXX` can override the compiler while preserving
the kernel ABI flags. Cache entries include the kernel build ID and compiler
flags, so a same-version kernel rebuild invalidates old compiled code.

For code calling an optional module API, install that module and explicitly
link its library, for example:

```sh
conda install -c CHANNEL plumed-mklib plumed-ext-ves
export PLUMED_MKLIB_LDFLAGS="$CONDA_PREFIX/lib/libplumed-ext-ves.so"
plumed mklib action.cpp
```

Use `.dylib` on macOS. The native dependency loads VES code; action visibility
still follows the loader rules described in the extension documentation.

## GPU validation and distribution handoff

[GPU validation notebook](gpu-validation.ipynb) installs final Linux x86-64 CI
artifacts on an NVIDIA runtime, including Colab. Upload the ordinary conda and
OpenACC artifact ZIPs from the same source revision, then run the notebook and
save its result logs. It exercises CUDA 12.9 and 13.0 and the OpenACC NVIDIA
target. It never overrides driver metadata: a driver too old for an output is
an incomplete matrix, not a successful GPU test. Alternatively index the artifact
packages as a local channel and run `run-gpu-test.sh` with `VERSION`,
`GPU_CONDA_CHANNEL` and `GPU_TEST_ROOT` set. The notebook needs a user-owned GPU
session; CPU CI does not establish GPU numerical correctness.

Generated recipes/catalogs are build products. Curate module requirements and
compatibility-package choices in their `conda.mk` files, regenerate with
`python extensions/conda_recipe.py`, and build the main package before its
outputs. Publish the exact base and dependent outputs together. The current
prototype is available as local/fork CI artifacts; public channel/feedstock
publication and the upstream merge are separate maintainer actions. OpenACC
release builders must provision the documented NVHPC SDK; users need no SDK.

Separate build jobs must reuse the exact kernel archive rather than rebuilding
another kernel with the same conda version/build label: native module build IDs
also capture the configured build. The focused CI shares one Linux kernel
artifact between CPU/CUDA and OpenACC jobs. After restoring and indexing that
archive, `PLUMED_CONDA_SKIP_BASE=yes` makes the build script use it; the script
requires the matching local archive and `CONDA_BLD_PATH`/`VERSION` settings.
