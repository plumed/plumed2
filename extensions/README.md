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
`PLUMED_PLUGIN_PATH` adds colon-separated registration directories (semicolons
on Windows), including
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

Action metadata commands (`manual`, `gentemplate --list`, `gen_json`, and
`gen_vim`) use the same runtime selection. For example,
`PLUMED_AUTOLOAD=0 plumed --load-module=ves gentemplate --list` lists the base
and VES actions, with their logical names. Installed Bash completion also
forwards explicit global `--load-module` selections into its inquiry, so commands
provided by optional modules can be completed. The installed JSON and Vim files
are build-time snapshots; use `gen_json` for a fresh JSON snapshot. Interactive
Vim support updates its own cache as described in `new-manual/vim.md`.

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
  --module cuda_coordination --variant cuda13.0 --features CUDA \
  --library ../../CudaCoordination.so --output plugin-registration
```

Compile `plugin-registration/module_identity.cpp` into the plugin library, link
to the selected base kernel and its native dependencies, then install the library
and `.plumed` file at the declared locations. Use the platform library suffix.
Separate conda recipes provide the Python, CUDA and OpenACC variants described
below. Ordinary CI validates compilation, loading and OpenACC host execution;
GPU numerical validation uses the channel-based Colab notebook after publication.

## Dependency inventory

The inventory and conda recipes, including `conda/plumed/meta.yaml`, are generated files and are not tracked. Run
`python3 extensions/conda_recipe.py` before invoking conda-build directly;
`conda/run-conda-build.sh` does this automatically. Configure and extension builds
regenerate the inventory themselves. The generators read `module.type`, exact
Makefile `USE` dependencies, source `__PLUMED_HAS_*` guards, module license files,
and module-local `conda.mk` declarations. See the complete
[configuration reference](../conda/module-options.md) for build tools, supported
platforms, variants, flags, feature policy and tests:

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
The base recipe is generated from `conda/plumed/meta.yaml.in` and the requirements
of `always` and `default-on` modules, including their native `USE` closure.
Platform selectors and explicit version/build constraints are retained; requirements
from default-off modules stay in their own outputs.
The extension recipe selects all default-off modules and opted-in empty compatibility
packages for modules supplied by the kernel. DRR declares Boost;
Fourier declares non-MPI FFTW; `isdb-pytorch` (EMMIVOX), PyTorch and metatomic declare CPU LibTorch and
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

`isdb` stays enabled in the default kernel without LibTorch. Its Torch-dependent
EMMIVOX action lives in the default-off `isdb-pytorch` module, installed with
`conda install plumed-ext-isdb-pytorch`. It can also be enabled in a source build
with `./configure --enable-modules=reset:+isdb-pytorch --enable-libtorch`, or built
separately with `--enable-extension-modules=isdb-pytorch`. Autoloading,
`LOAD MODULE=isdb-pytorch` and per-module feature queries follow the same
contract as other extensions. It does not depend on the PLUMED `pytorch` module.

The three existing plugins use the common registration mechanism and have
local conda declarations and generated recipes, described below.

## Installed numerical checks

`conda/plumed-extensions/regtests.json` maintains the required numerical cases
and the extension packages that supply their native dependencies. The conda
harness installs those packages at the base version before running the cases.
The current selection covers VES static bias/forces, Fourier grids, DRR
checkpoint/restart, PyTorch strided tensors/derivatives and EMMIVOX
fixed-B-factor forces/derivatives and B-factor checkpoint/restart.
The full EMMIVOX fixture is retained for native tests; its sampled B-factor
derivatives can differ across LibTorch backends because of single-precision
cancellation near zero density residuals. The conda selection keeps strict
reference checks for the stable calculation and bitwise restart checks for sampling. Symfunc is declared
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
3.11–3.15 variants declared in `plugins/pycv/conda.mk`. Every variant
links to and pins the same main `plumed` package, without bundling another kernel.
Python and NumPy are plugin runtime requirements; pybind11 is needed only to build.

Like other extensions, pycv loads automatically. With `PLUMED_AUTOLOAD=0`, use
`LOAD MODULE=pycv` in the input or `plumed --load-module=pycv ...` on the command
line. `plumed modules --has pycv` checks runtime availability. `python -m pycv` continues to report the installed library path. Use
`LOAD MODULE=pycv` with autoloading disabled so its private runtime dependencies load with it.

`extensions/conda_recipe.py` also generates `conda/plumed-pycv/meta.yaml`, using
its module-local curated declarations and license. The build installs both
Python components together and registers the library under `lib/plumed/plugins`.
On macOS the same package also installs a private `PythonCVRuntime.dylib`
dependency, so Python-hosted and standalone runs share the right interpreter
symbols. It is not registered as a module: `plumed modules --list` shows only
`pycv`, and users never need a second `LOAD MODULE` command. No extra kernel or compiler is required at runtime.
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

### OpenACC parallel task manager plugin

`plumed-ext-openacc-ptm` is built against the exact main kernel with an external
NVHPC SDK. The initial recipe targets Linux x86-64 with NVHPC 25.11/CUDA 13.0.
It produces one library containing NVIDIA GPU and serial-host targets, with
machine code for multiple GPU architectures. These are runtime choices within
one package, not alternative device-specific packages. AMD and Apple GPUs are
not supported by this compiler backend.

The native module name is `openacc_ptm`. It autoloads normally, or can be loaded
with `LOAD MODULE=openacc_ptm` / `plumed --load-module=openacc_ptm ...` when
`PLUMED_AUTOLOAD=0`. `plumed modules --has openacc_ptm` checks availability and
`plumed modules --has-feature openacc` checks the visible compiled feature.

Select the execution target with `NVCOMPILER_ACC_DEVICE_TYPE=HOST` or `NVIDIA`.
The recipe statically links the NVIDIA runtime, includes its license notices,
and adds no SDK dependency to the main kernel. A compatible NVIDIA driver is
needed only for GPU execution; the package can install and run host tests on a
machine without one. Multicore execution is not included because its NVIDIA
OpenMP runtime needs separate compatibility validation with the GCC-built base.

To build, provide an isolated SDK installation, generate the recipes and build
against a locally available main package:

```sh
bash conda/fetch-nvhpc.sh 25.11 /path/to/private-sdk
export NVHPC_HOME=/path/to/private-sdk/opt/nvidia/hpc_sdk
export VERSION=2.11.0dev0
python extensions/conda_recipe.py
conda-build -c local -c conda-forge conda/plumed-openacc
```

The fetch helper downloads the official SDK, verifies its published SHA256 and
extracts compiler/CUDA components without root or system installation changes.
An existing matching SDK can be used instead. The compiler is a build tool only.
The ordinary conda build script also builds this recipe when `NVHPC_HOME` is set.

The package test runs existing DISTANCE value/force and SAXS/SANS chi-squared derivative references on the installed
multi-target binary with host execution, testing both autoload and explicit
loading. To require GPU execution of those same comparisons, run:

```sh
python extensions/tests/openacc.py --prefix "$CONDA_PREFIX" --device nvidia
```

GPU mode requires `nvidia-smi` and checks NVIDIA runtime kernel-launch messages,
so CPU fallback cannot count as a GPU test. Host results do not validate device
code or transfers. `bash conda/run-openacc-test.sh` separately checks fresh
installation/removal, absence of runtime compiler dependencies, and unchanged
main kernel bytes.

### Modules moved into the kernel

Set `CONDA_DUMMY_PACKAGE = yes` in a kernel-built module's `conda.mk` to retain
its historical `plumed-ext-<name>` package as an empty compatibility output.
The generated package pins the exact main kernel and installs no module library
or registration. Its test requires the module in the static kernel config.
OPES uses this setting so existing `plumed-ext-opes` installations can upgrade
when OPES is supplied by the main package. No shortcut dependencies are added.

### Per-module build features and regression prerequisites

`plumed modules --module drr --has-feature boost_serialization -q` tests a
feature in DRR's build. Kernel-built modules report the kernel feature vector;
separate modules report their own final enabled compile-time features. Features
from other modules or dependencies are not unioned into a scoped query. Without
`--module`, the existing aggregate inquiry remains available. `--module NAME
--list` lists that visible module's enabled features. Queries do not activate a
module: use normal autoloading or pre-query `--load-module=NAME` when needed.

Regression configs declare `plumed_needs="drr:boost_serialization"`.
The runner
uses scoped inquiries when supported, otherwise checks both static module
presence and the feature through `plumed config`. A false modern inquiry never
falls back to kernel-wide flags. Requirements are checked before reading the
input, so a `LOAD MODULE` only inside that input cannot satisfy a prerequisite.
The canonical feature is `boost_serialization`; there is no `boost` alias.

### Optional runtime C++ tools

`conda install plumed-mklib` adds the runtime compiler toolchain, `plumed mklib`
and `LOAD FILE=MyAction.cpp` support to the exact main package. The base keeps
its development headers/build settings but does not ship the mklib script or
launcher. The tooling package adds no PLUMED module. Activate the environment
to select its compiler/SDK, or use `PLUMED_MKLIB_CXX` for an explicit compatible
compiler. Additional native module/library links use `PLUMED_MKLIB_LDFLAGS`.

Generated libraries link the kernel and compiler runtimes; imported inline
BLAS/LAPACK calls additionally link their actual provider. Other private kernel
libraries do not require development linker files just to compile a basic
action. C++ LOAD cache filenames include the kernel build ID; direct mklib -n
also records kernel identity, compiler/link flags and the compiler helper content. Existing atomic publication
and simultaneous-builder handling are retained.
