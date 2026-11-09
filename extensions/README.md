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
The build also inspects undefined object symbols for BLAS/LAPACK calls emitted
by public inline matrix helpers. It links the matching numerical provider flags
from the base and records their `__PLUMED_HAS_EXTERNAL_*` macros. Each installed
module owns an accompanying JSON dependency report beside its registration;
`_extensions/native-dependencies.json` contains the complete build report.

Native package mappings are based on the conda-forge feedstocks for
[Boost](https://github.com/conda-forge/boost-feedstock),
[FFTW](https://github.com/conda-forge/fftw-feedstock),
[LibTorch](https://github.com/conda-forge/pytorch-cpu-feedstock),
[metatomic](https://github.com/conda-forge/libmetatomic-torch-feedstock), and
[HDF5](https://github.com/conda-forge/hdf5-feedstock), with the earlier source
inventory audited on 2026-11-08. They are availability mappings, not evidence
that extension builds with these dependencies have passed.

## Initial conda coverage

The initial recipe provides 30 optional modules: the 29 with no mandatory
module-specific native feature plus Fourier. The existing conda base already
provides non-MPI FFTW, so Fourier adds no new runtime stack. Its extension links
FFTW explicitly instead of importing unrelated private kernel libraries.

| Module or plugin | Required feature macros | Native dependency | Conda-forge mapping | Delivery |
| --- | --- | --- | --- | --- |
| fourier | `__PLUMED_HAS_FFTW` | `fftw` (non-MPI) | Available | Initial recipe |
| drr | `__PLUMED_HAS_BOOST_SERIALIZATION` | `libboost-devel` / `libboost` | Available | Later native-dependency builds |
| pytorch | `__PLUMED_HAS_LIBTORCH` | `libtorch` | Available | Later CPU/CUDA variants |
| metatomic | `__PLUMED_HAS_LIBMETATOMIC`, `__PLUMED_HAS_LIBTORCH` | Native metatomic/metatensor Torch libraries and LibTorch | Available | Later native-dependency builds |
| ttsketch | `__PLUMED_HAS_LIBITENSOR`, `__PLUMED_HAS_LIBHDF5` | ITensor v3 C++ and HDF5 | HDF5 available; ITensor recipe not verified | Deferred |
| kenref | `__PLUMED_HAS_KENREF` | `kenref_core` and interface headers | Recipe not verified | Deferred |
| cuda_coordination plugin | `__PLUMED_HAS_CUDA` | CUDA compiler/runtime | Available on supported CUDA platforms | Later variants |
| openacc_ptm plugin | `__PLUMED_HAS_OPENACC` | NVIDIA HPC SDK compiler/runtime | Compiler recipe not verified | Deferred |
| pycv plugin | `__PLUMED_HAS_PYTHON_EMBED` | Python development API and pybind11 | Available | Later plugin recipe |

Availability denotes package/recipe mapping, not a validated platform/version
matrix. All plugins can use the common identity/registration mechanism now.
The initial conda recipe does not yet package the three existing plugin builds.
