Platforms and features
======================

PLUMED can be installed in three ways:

* a **standard build** from the sources, compiled with `./configure` and `make`
  on Linux and macOS (the official build system), or with CMake and Visual Studio
  on Windows (experimental);
* the **Python package** (`pip install plumed`), whose wheels bundle a PLUMED
  kernel and the `plumed` command;
* the **conda packages**.

The tables below summarize which features are available in each case. A feature
marked as "no" is not compiled in, and the related actions or command line tools
report an error when they are used. The features of an installed PLUMED are
listed by `plumed config show`. How to install PLUMED is explained in
[README.md](README.md#required-software) and in the
[user documentation](https://www.plumed.org/doc).

Standard builds
---------------

| | Linux | macOS | Windows |
|---|---|---|---|
| Build system | `./configure` and `make` | `./configure` and `make` | CMake and Visual Studio (experimental) |
| Libraries | shared kernel, static archive, wrapper library | same | static kernel (to be linked into an MD code), or shared kernel with `-DPLUMED_SHARED_KERNEL=ON` |
| Modules (see [README.md](README.md#modules-and-plugins)) | any (`--enable-modules`) | any | `always` and `default-on` ones; VES with `-DPLUMED_ENABLE_VES=ON` |
| MPI | external MPI library, or embedded thread-MPI | same | embedded thread-MPI only |
| OpenMP | yes | yes (with Apple clang, needs libomp) | no |
| Expressions compiled with asmjit | yes (x86_64 and 64 bit ARM) | yes | yes |
| Optional libraries (BLAS/LAPACK, zlib, HDF5, FFTW, chemfiles, LibTorch, metatomic, ...) | yes, if found or requested | same | HDF5 and chemfiles (`-DPLUMED_ENABLE_HDF5=ON`, `-DPLUMED_ENABLE_CHEMFILES=ON`); internal BLAS and LAPACK |
| Plugins and your own actions (`LOAD`, `plumed --load`, `plumed mklib`), `PLUMED` action with `KERNEL` | yes | yes | no |
| Shared kernel loaded at run time (`PLUMED_KERNEL`, Python package) | yes | yes | yes, with `-DPLUMED_SHARED_KERNEL=ON` |
| Patching MD codes (`plumed patch`) | yes | yes | no |
| Tools written as shell scripts | yes | yes | no |
| Running external programs | yes | yes | no |
| Python package compiled from the sources | yes | yes | yes |
| Stack trace in error messages | yes | yes | yes (function names need the debug information, e.g. `-DPLUMED_MSVC_DEBUG_INFO=ON`) |
| File names with non-ASCII characters | yes | yes | yes (passed to the kernel as UTF-8; some contributed modules only support the system code page) |
| Regression tests (`make test`) | all | all | a subset, run with a script (`.ci/windows-regtests`) |


Python wheels
-------------

All the wheels are built from the same sources with the same choice of modules
and libraries (HDF5 and chemfiles are linked statically and their symbols are
not exported, so that they do not interfere with other copies of them loaded in
the same process, e.g. by h5py or by an MD code), so that a simulation input gives the same results on every
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
| HDF5 files (`.h5`), more trajectory formats with chemfiles (driver `--itraj`, GROMACS `.tpr` files) | yes | yes | yes |
| Compressed files (zlib), FFTW | no | no | no |
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

Conda packages
--------------

The recipes in the [conda](conda/README.md) directory build these packages
from the same sources, for the `plumed` channel:

* `plumed`: the kernel, the wrappers, the headers and the `plumed` command;
* `plumed-mpi` (optional): an MPICH or Open MPI provider for the same kernel,
  and the `plumed-mpi` executable;
* `py-plumed`: the Python interface, using the kernel of `plumed`.

| | Linux x86_64, aarch64 | macOS x86_64, arm64 | Windows |
|---|---|---|---|
| Packages | `plumed`, `plumed-mpi`, `py-plumed` | same | none (use the Python wheels or a standard build) |
| Kernel built with | `./configure` and `make`, conda compilers | same | |
| Libraries | shared kernel, wrapper library (no static archive) | same | |
| Modules (see [README.md](README.md#modules-and-plugins)) | `always` and `default-on` ones | same | |
| MPI | embedded thread-MPI; MPICH or Open MPI with `plumed-mpi`, without rebuilding the kernel | same | |
| OpenMP | yes | yes (llvm-openmp) | |
| Expressions compiled with asmjit | yes | yes | |
| BLAS and LAPACK | from conda (`libblas`, `liblapack`) | same | |
| Compressed files (zlib), FFTW, HDF5 files | yes | yes | |
| chemfiles (`plumed driver --itraj`) | yes (no GROMACS TPR files on aarch64, where only chemfiles 0.10.3 is available) | yes | |
| Python interface | `py-plumed`, CPython 3.11 and later | same | |
| Plugins and your own actions (`LOAD`, `plumed mklib`; compiling them needs a C++ compiler) | yes | yes | |
| Patching MD codes (`plumed patch`) and tools written as shell scripts | yes | yes | |
| Linking C/C++ programs and MD codes | yes, as with a standard build | same | |

The packages on [conda-forge](https://anaconda.org/conda-forge/plumed)
(`conda install -c conda-forge plumed`) are built from a separate recipe, and
their features can differ: at the moment, each MPI implementation (`nompi`,
`mpich`, `openmpi`) has its own build of the whole kernel, and HDF5 and
chemfiles are not included.
