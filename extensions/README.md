# Optional native modules

One `plumed` base supports separately installed `plumed-ext-<module>` packages.
Each extension is rebuilt for every base update and pinned to its exact version
and conda build. Serial, MPICH and Open MPI use the same extension through mpiw.
Only hard native dependencies are honored; shortcut dependencies are not tracked.

## Build against an installed base

Use the same source version and a compatible compiler/toolchain as the base:

```sh
./configure --with-extension-base=/path/to/base \
  --enable-extension-modules=ves:pamm --prefix=/path/to/extensions
make
make install
```

This builds VES, PAMM and its hard adjmat dependency; it does not rebuild or
install a kernel. Python 3 is needed for this build interface, not at runtime.
The installed base supplies its compiler flags and public headers. An explicitly
set CXX/CXXFLAGS/CPPFLAGS can adapt paths/toolchain activation for a package build;
it remains the builder's responsibility to preserve C++ ABI compatibility.
Available optional native features are detected with compile/link probes against
this installation, independently of the base's feature macros. For this loading
proof of concept every selected module is built even when an optional library
is unavailable. Such a module can offer fewer actions; developers will later
classify which libraries are compulsory.

Each package owns its library in `PREFIX/lib` and a small `.plumed` registration
file in `PREFIX/lib/plumed/plugins`. The file declares `module`, exact `build`,
`library` (relative to the registration directory), optional `variant`, and
`requires` (direct hard module dependencies). Compatibility of the complete
selected dependency closure is checked before loading it. Its native shared
library links must independently resolve every dependency.

## Loading

The directory beside the actual loaded kernel is discovered automatically.
`PLUMED_PLUGIN_PATH` adds colon-separated registration directories, including
separate HPC installations. HPC modulefiles must expose the registration
directories of hard module dependencies too, typically through modulefile
dependencies; native library search paths alone do not supply registration metadata.
Empty elements do not search the current directory.
Discovery is cached for this kernel; install/remove changes take effect in a new
process. A missing directory is normal. Broken selected registrations or failed
loads are errors. Recursive extension activation from library initialization
is rejected with the name of the module being loaded.

Disable eager loading with `PLUMED_AUTOLOAD=0`; named loading still works:

```plumed
LOAD MODULE=ves
```

```sh
PLUMED_AUTOLOAD=0 plumed --load-module=ves driver --plumed plumed.dat --noatoms
plumed modules --has ves
plumed modules --list
```

Explicit loading affects only the requesting object or CLI context. Hard
dependencies become available silently in that context. Repeated loads, including
already-built-in modules, succeed without changing visibility elsewhere.
Autoloaded modules are shared by the kernel and retained across object lifetimes.
`plumed modules --list` reports built-in and successfully loaded modules.
`plumed modules --has ves pamm` checks all requested modules and returns 0 only
when all are visible; `-q` suppresses output. Runtime native feature checks use
`plumed modules --has-feature libtorch`; they combine the base features with
features exported by visible extension libraries. This reports availability in
this invocation, including silently loaded dependencies, rather than installation
alone. `PLUMED_AUTOLOAD=0` affects these queries normally.

Both `plumed config` and standalone `plumed-config` describe the static base
build, including module and library features. The standalone script can be used
when cross compiling without executing the target kernel. For example, installing
PyTorch can make `plumed modules --has-feature libtorch` succeed while
`plumed config has libtorch` still describes a base built without Torch.

Two variants with the same module identity are an error. HPC modulefiles should
expose only the selected variant's directory. For troubleshooting,
`PLUMED_MODULE_OVERRIDE='name=/absolute/library;other=/absolute/library'` selects
exact libraries among registrations; it does not bypass compatibility checks or
make mutually incompatible variants usable together. Native dependency search
paths remain the responsibility of the library build and HPC modulefile.

## Existing plugins

Plugins use the same registration format and loading contract. Generate an
identity translation unit and registration file before the plugin's own build:

```sh
python3 extensions/register.py --plumed /path/to/base/bin/plumed \
  --module cuda_coordination --variant cuda13 \
  --library ../../libplumed-ext-cuda-coordination.so --output plugin-registration
```

Compile `plugin-registration/module_identity.cpp` into the plugin library, link
to the selected base kernel and its native dependencies, then install the library
and `.plumed` file at the declared locations. Use the platform library suffix.
CUDA/OpenACC numerical and variant validation requires suitable hardware/tools;
these conda variants are a subsequent stage, not supplied by the initial recipe.

## Dependency inventory

The inventory and conda recipe are generated files and are not tracked. Run
`python3 extensions/conda_recipe.py` before invoking conda-build directly;
`conda/run-conda-build.sh` does this automatically. Configure and extension builds
regenerate the inventory themselves. The generators read `module.type`, exact
Makefile `USE` dependencies, source `__PLUMED_HAS_*` guards, module license files,
and module-local `conda.mk` declarations:

```make
CONDA_HOST_REQUIREMENTS = libboost-devel
CONDA_RUN_REQUIREMENTS = libboost
# Libraries without a verified conda-forge package are recorded, not required:
CONDA_UNAVAILABLE_LIBRARIES = example-library
```

Quote MatchSpecs containing spaces. Conda selectors are supported at the end of
a declaration, for example `# [linux]`. These files are manually curated inputs and are never overwritten by the
generators. They declare package policy;
compiler probes only decide which feature macros and native links to enable.
The current recipe selects all 35 default-off modules. DRR declares Boost;
Fourier declares non-MPI FFTW; PyTorch and metatomic declare CPU LibTorch and
native metatomic/metatensor libraries; TTsketch builds against non-MPI HDF5 but permits ABI-compatible MPI variants
at runtime, while its
unavailable ITensor dependency does not exclude it. KEnRef is also retained.
Each package gets the license and notices from its own module directory, with
the root license used only when a module has no notice.

The build inspects undefined object symbols for BLAS/LAPACK calls emitted by
public inline matrix helpers and links the matching numerical provider from
the base. Installed JSON reports beside each registration list enabled optional
features, native libraries and direct hard module dependencies. The full report
is `_extensions/native-dependencies.json`. No shortcut dependencies are inferred.

The three existing plugins use the common registration mechanism and also have
local conda declarations. Plugin-specific recipes and CUDA variants remain a
subsequent stage.

## Installed numerical checks

`conda/plumed-extensions/regtests.json` maintains the required numerical cases
and the extension packages that supply their native dependencies. The conda
harness installs those packages at the base version before running the cases.
The current selection covers VES static bias/forces, Fourier grids, DRR
checkpoint/restart and PyTorch strided tensors/derivatives. Symfunc is declared
only for the Fourier test input, without becoming a Fourier package dependency.

The selection can explicitly assert its prerequisites with `force_requirements`.
The runner sets `PLUMED_REGTEST_FORCE_REQUIREMENTS=yes` only for those cases,
bypassing module/library eligibility checks in the ordinary regression harness.
Actual calculations, hooks and reference comparisons still run unchanged; a
missing action, numerical mismatch or remaining skip fails the package check.
Ordinary regression runs keep their normal eligibility checks. This curated
validation does not depend on choosing a per-module feature-query syntax.

Plugin builders can pass `--features LIBRARY_MACRO` to `extensions/register.py`.
The identity translation unit exports those `__PLUMED_HAS_LIBRARY_MACRO` flags
only when enabled by its compile flags, following the same runtime inquiry
contract as source modules.

### Python collective-variable plugin

`conda install plumed-ext-pycv` installs the pycv plugin and its companion Python
extension for the environment's Python minor version. The recipe builds Python
3.11–3.15 variants from `conda/plumed-pycv/conda_build_config.yaml`. Every variant
links to and pins the same main `plumed` package, without bundling another kernel.
Python and NumPy are plugin runtime requirements; pybind11 is needed only to build.

Like other extensions, pycv loads automatically. With `PLUMED_AUTOLOAD=0`, use
`LOAD MODULE=pycv` in the input or `plumed --load-module=pycv ...` on the command
line. `plumed modules --has pycv` checks runtime availability. `python -m pycv` continues to report the installed library path. Use
`LOAD MODULE=pycv` with autoloading disabled so its runtime helpers load first.

`extensions/conda_recipe.py` also generates `conda/plumed-pycv/meta.yaml`, using
its module-local curated declarations and license. The build installs both
Python components together and registers the library under `lib/plumed/plugins`.
On macOS the same package also installs a `pycv_runtime` helper module, loaded
as a dependency before pycv, so Python-hosted and standalone runs share the
right interpreter symbols. No extra kernel or compiler is required at runtime.
Its package tests run the existing Python-hosted suite and numerical driver
fixtures for both autoloading and explicit module loading.

### CUDA coordination plugin

`plumed-ext-cuda-coordination` has CUDA 12.9 and 13.0 builds on Linux. Select one
with `cuda-version=12.9` or `cuda-version=13.0`; they are variants of the same
package, so conda installs one at a time. Each links to the exact main `plumed`
version/build and the shared CUDA runtime. The kernel has no CUDA compiler or
runtime dependency. Compiler tools and CCCL headers are build requirements only.

The plugin autoloads as `cuda_coordination`. With autoloading disabled use
`LOAD MODULE=cuda_coordination` or `plumed --load-module=cuda_coordination ...`.
The loader rejects multiple variants exposed through `PLUMED_PLUGIN_PATH`,
including variants outside the active conda prefix.

Package tests without a GPU verify native loading, both action registrations,
feature exports and explicit loading. The install/remove harness switches CUDA
variants and checks that the main kernel bytes stay unchanged. Driver metadata
is overridden only for those CPU-only CI checks. Actual GPU numerical validation
can be requested with `python extensions/tests/cuda_coordination.py --prefix
"$CONDA_PREFIX" --cuda-version 12.9 --gpu`; it runs existing double/float value
and derivative fixtures in autoload and explicit module modes and rejects skips.
Passing loading checks alone does not validate GPU calculations.
