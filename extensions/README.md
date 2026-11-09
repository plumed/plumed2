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
plumed config module ves
plumed --list-modules
```

Explicit loading affects only the requesting object or CLI context. Hard
dependencies become available silently in that context. Repeated loads, including
already-built-in modules, succeed without changing visibility elsewhere.
Autoloaded modules are shared by the kernel and retained across object lifetimes.
`config module` reports built-in and successfully loaded modules. Other config
fields, including `makefile_conf`, describe the base build. The standalone
`plumed-config` script remains a static base-build inquiry for cross compilation;
use the native `plumed config module` command for runtime availability.

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
a declaration, for example `# [linux]`. These files declare package policy;
compiler probes only decide which feature macros and native links to enable.
The current recipe selects all 35 default-off modules. DRR declares Boost;
Fourier declares non-MPI FFTW; PyTorch and metatomic declare CPU LibTorch and
native metatomic/metatensor libraries; TTsketch declares non-MPI HDF5 while its
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
