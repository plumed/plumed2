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
Modules needing additional external-library features are rejected when the base
configuration does not provide them, rather than producing nonfunctional stubs.

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

`python3 extensions/inventory.py` regenerates `catalog.json` from `module.type`,
Makefile `USE` and source preprocessor feature guards. It records
`__PLUMED_HAS_*` evidence separately from reviewed mandatory native dependencies.
Optional OpenACC paths do not make those whole modules require a GPU.
Conda-forge recipe availability does not establish a compatible version/platform
solve. ITensor, KEnRef and the NVIDIA OpenACC compiler need further packaging
verification. No shortcut dependency inventory is maintained.

Native package mappings are based on the conda-forge feedstocks for
[Boost](https://github.com/conda-forge/boost-feedstock),
[FFTW](https://github.com/conda-forge/fftw-feedstock),
[LibTorch](https://github.com/conda-forge/pytorch-cpu-feedstock),
[metatomic](https://github.com/conda-forge/libmetatomic-torch-feedstock), and
[HDF5](https://github.com/conda-forge/hdf5-feedstock), with the earlier source
inventory audited on 2026-11-08. They are availability mappings, not evidence
that extension builds with these dependencies have passed.
