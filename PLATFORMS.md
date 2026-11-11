Platforms and features
======================

PLUMED can be installed in two ways:

* a **standard build** from the sources, compiled with `./configure` and `make`
  on Linux and macOS (the official build system), or with CMake and Visual Studio
  on Windows (experimental);
* the **Python package**, `python -m pip install plumed`, whose wheels bundle a
  PLUMED kernel and the `plumed` command.

The tables below summarize which features are available in each case. A feature
marked as "no" is not compiled in, and the related actions or command line tools
report an error when they are used.

Standard builds
---------------

| | Linux | macOS | Windows |
|---|---|---|---|
| Build system | `./configure` and `make` | `./configure` and `make` | CMake and Visual Studio (experimental, see `CMakeLists.txt`) |
| Libraries | shared kernel, static archive, wrapper library | same | static kernel (to be linked into an MD code), or shared kernel with `-DPLUMED_SHARED_KERNEL=ON` |
| Modules (see [README.md](README.md#modules-and-plugins)) | any (`--enable-modules`) | any | `always` and `default-on` ones; VES with `-DPLUMED_ENABLE_VES=ON` |
| MPI | external MPI library, or embedded thread-MPI | same | embedded thread-MPI only |
| OpenMP | yes | yes (with Apple clang, needs libomp) | no |
| Expressions compiled with asmjit | yes (x86_64 and 64 bit ARM) | yes | yes |
| Optional libraries (BLAS/LAPACK, zlib, HDF5, FFTW, chemfiles, LibTorch, metatomic, ...) | yes, if found or requested | same | no (internal BLAS and LAPACK only) |
| Plugins and your own actions (`LOAD`, `plumed --load`, `plumed mklib`), `PLUMED` action with `KERNEL` | yes | yes | no |
| Shared kernel loaded at run time (`PLUMED_KERNEL`, Python package) | yes | yes | yes, with `-DPLUMED_SHARED_KERNEL=ON` |
| Patching MD codes (`plumed patch`) | yes | yes | no |
| Tools written as shell scripts | yes | yes | no |
| Running external programs | yes | yes | no |
| Python package compiled from the sources | yes | yes | yes |
| Stack trace in error messages | yes | yes | yes (function names need the debug information, e.g. `-DPLUMED_MSVC_DEBUG_INFO=ON`) |
| File names with non-ASCII characters | yes | yes | yes (passed to the kernel as UTF-8; some contributed modules only support the system code page) |
| Regression tests (`make test`) | all | all | a subset, run with a script (`.ci/windows-regtests`) |

Building on Windows:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release --parallel
    cmake --install build --prefix <prefix>


Python wheels
-------------

All the wheels are built from the same sources with the same choice of modules
and libraries, so that a simulation input gives the same results on every
platform. They are meant to analyze data from Python and to run small
simulations; for production runs coupled to an MD code, use a standard build
(or a conda package) configured for your machine.

| | Linux x86_64, aarch64 | macOS x86_64, arm64 | Windows x86_64 |
|---|---|---|---|
| Minimum OS | glibc 2.17 (manylinux2014) | macOS 10.15 (Intel), 11.0 (Apple Silicon) | 64 bit Windows |
| Python | CPython 3.11 and later | CPython 3.11 and later | CPython 3.11 and later |
| Kernel built with | `./configure` and `make`, GCC | `./configure` and `make`, Apple clang | CMake, Visual Studio (static C/C++ runtime) |
| Modules (see [README.md](README.md#modules-and-plugins)) | `always` and `default-on` ones | same | same |
| Expressions compiled with asmjit | yes | yes | yes |
| Parallelism | embedded thread-MPI (`--tmpi-threads N`); no MPI library, no OpenMP | same | same |
| BLAS and LAPACK | internal | internal | internal |
| Compressed files (zlib), FFTW, HDF5 files, chemfiles | no | no | no |
| Python interface (`Plumed`, `InputBuilder`, I/O helpers) | yes | yes | yes |
| `plumed` command line tools written in C++ (`driver`, `sum_hills`, ...) | yes | yes | yes |
| Tools written as shell scripts (`plumed patch`, `config`, `mklib`, `newcv`, ...) | yes | yes | no |
| Plugins and your own actions (`LOAD`, `plumed mklib`; compiling them needs a C++ compiler) | yes | yes | no |
| Running external programs (e.g. `MOLINFO` selections with `PYTHON_BIN`) | yes | yes | no |
| Stack trace in error messages | yes | yes | without function names (no debug information) |
| File names with non-ASCII characters | yes | yes | yes |
| Linking C/C++ programs to the bundled kernel | yes (headers, `compile_options.sh`, pkg-config) | same | C interface only (`plumedKernel.lib`, pkg-config) |

Windows on ARM, 32 bit systems and Linux distributions based on musl (e.g.
Alpine) have no wheels: there, `pip` builds the package from the sources.
The Python package can also be installed without its kernel and used with a
kernel from a standard build (see [python/README.rst](python/README.rst)).

