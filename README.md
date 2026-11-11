[![Homepage](https://img.shields.io/badge/Home-plumed.org-green.svg)](https://www.plumed.org)
[![Google group](https://img.shields.io/badge/Google_group-plumed--users-green.svg)](https://groups.google.com/g/plumed-users)
[![codecov](https://codecov.io/gh/plumed/plumed2/branch/master/graph/badge.svg)](https://codecov.io/gh/plumed/plumed2)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Github Releases](https://img.shields.io/github/release/plumed/plumed2.svg)](https://github.com/plumed/plumed2/releases)
[![MacPorts package](https://repology.org/badge/version-for-repo/macports/plumed.svg)](https://repology.org/project/plumed/versions)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/plumed/badges/version.svg)](https://anaconda.org/conda-forge/plumed)
[![AUR package](https://repology.org/badge/version-for-repo/aur/plumed.svg)](https://repology.org/project/plumed/versions)
[![DPorts package](https://repology.org/badge/version-for-repo/dports/plumed.svg)](https://repology.org/project/plumed/versions)
[![FreeBSD port](https://repology.org/badge/version-for-repo/freebsd/plumed.svg)](https://repology.org/project/plumed/versions)
[![Spack package](https://repology.org/badge/version-for-repo/spack/plumed.svg)](https://repology.org/project/plumed/versions)
[![Bluesky](https://img.shields.io/bluesky/followers/plumed.org?label=Follow)](https://bsky.app/profile/plumed.org)

Branches and releases
---------------------

Several branches and tags are stored on the git repository.

Branches named `v2.X` correspond to release branches.

Master branch may contain non tested features and is not expected to be used by non-developers.
It typically contains features that will be available on the next release.

Tags named `v2.XbY` correspond to beta releases, use it with care.
Tags named `v2.X.Y` correspond to official releases, use the latest available.

In addition, the repository contains a number of other branches related to specific features.
Please contact the developers that are committing on those branches before basing your work
there, since they might contain temporary work and might be rebased later.
For instance, branch `testdoc` is setup so as to push a test copy of the manual
and is often force pushed.

To report problems found on beta or official releases, use the normal
[plumed-users@googlegroups.com](mailto:plumed-users@googlegroups.com)
mailing list. Please state exactly which version you are using.
To report problems found on `master` branch, use the
[plumed2-git@googlegroups.com](mailto:plumed2-git@googlegroups.com) mailing list.
This is also the correct place for discussions about new features etc.
When reporting please provide the git hash (you can obtain it with `git rev-parse HEAD`).

Status
------

Below you find the status on [GitHub Actions](https://github.com/plumed/plumed2/actions) for the release branches.

| Branch   |      Status   | First stable release (year) | Still supported |
|:--------:|:-------------:|:--------:|:------:|
| master   | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/plumed/plumed2/actions/workflows/ci.yml) | 2026 (expected) | / |
| v2.10    | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=v2.10)](https://github.com/plumed/plumed2/actions/workflows/ci.yml)| 2025 | yes |
| v2.9     | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=v2.9)](https://github.com/plumed/plumed2/actions/workflows/ci.yml)   | 2023 | yes |
| v2.8     | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=v2.8)](https://github.com/plumed/plumed2/actions/workflows/ci.yml) | 2022 | no |
| v2.7     | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=v2.7)](https://github.com/plumed/plumed2/actions/workflows/ci.yml) | 2020 | no |
| v2.6     | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=v2.6)](https://github.com/plumed/plumed2/actions/workflows/ci.yml) | 2019 | no |
| v2.5     | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=v2.5)](https://github.com/plumed/plumed2/actions/workflows/ci.yml) | 2018 | no |
| v2.4     | [![CI](https://github.com/plumed/plumed2/actions/workflows/ci.yml/badge.svg?branch=v2.4)](https://github.com/plumed/plumed2/actions/workflows/ci.yml) | 2017 | no |
| v2.3     | Not available | 2016 | no |
| v2.2     | Not available | 2015 | no |
| v2.1     | Not available | 2014 | no |
| v2.0     | Not available | 2013 | no |

Content
-------

Here's a description of the content of each file and directory in the root PLUMED directory.

    .ci              : helper scripts used by continuous integration
    CHANGES          : change log
    CMakeLists.txt   : experimental CMake build for MSVC/Windows (not the official build system)
    COPYING.LESSER   : license
    Makefile         : makefile
    Makefile.conf.in : template configuration makefile
    PEOPLE           : list of authors
    PLATFORMS.md     : features available on each platform and in the Python package
    README.md        : this file
    VERSION.txt      : version file
    astyle           : a local version of astyle, used to format code
    cmake-msvc       : support files for the experimental CMake build
    conda            : conda packaging scripts
    configure        : configuration script
    configure.ac     : configuration script (autoconf)
    developer-doc    : developer documentation
    docker           : directory where Docker is generated
    fortran          : Fortran interfaces
    json             : tools to generate the JSON syntax file
    macports         : directory where Portfiles are generated
    new-manual       : user manual (MkDocs)
    patches          : patch scripts
    plugins          : optional plugins (e.g. pycv)
    python           : python package (wrappers and pip packaging)
    pythonformatter.sh : formats the python code (as astyle for C++)
    regtest          : regression tests, including reference results
    release.sh       : developer utility to publish releases
    ruff.toml        : configuration for the ruff Python linter
    scripts          : shell tools
    sourceme.sh.in   : template configuration script
    src              : source code
    stamp-h.in       : timestamp file used by autoconf
    test             : examples
    user-doc         : old tutorials and pages from the Doxygen manual, not part of the current manual
    vim              : directory where vim syntax is generated

Required software
-----------------

The following applies to the standard build on Linux and macOS. The features available 
on different architectures are described in: [PLATFORMS.md](PLATFORMS.md).

Required software:

* GNU make.
* C/c++ compiler (c++17 support is required as of version 2.10).
* A modern version of the `patch` command line tool.

Suggested software (libraries are checked by `./configure` and enabled if available):

* An external MPI library for distributed-memory parallel simulations, matching the MPI implementation used by your MD code. Serial builds can instead use embedded thread-MPI for local ranks; see [MPI runtime selection](src/mpiwrap/README.md).
* Optimized blas and lapack libraries (e.g. OpenBLAS or MKL). They are automatically replaced by an internal version if not available. They only make a difference for actions working on large matrices (e.g. DIAGONALIZE, INVERT_MATRIX or dimensionality reduction): the internal version is as fast as an optimized library for small matrices and faster than the reference netlib implementation.
* [Zlib library](https://zlib.net/) to use compressed data files.
* [HDF5 library](https://www.hdfgroup.org/solutions/hdf5/) to read and write tables in HDF5 files (`.h5`).
* [FFTW library](https://www.fftw.org/) for FOURIER_TRANSFORM (`fourier` module).
* [Chemfiles library](https://chemfiles.org/) to read more trajectory formats in `plumed driver` (only used if requested with `--enable-libchemfiles`).
* [MkDocs](https://www.mkdocs.org) to build the user manual (Python requirements are listed in `new-manual/requirements.txt`).
* [Doxygen](https://www.doxygen.nl) to build the developer manual. Doxygen might need:
  * [Graphviz](https://www.graphviz.org) to show class hierarchy in
    developer manual.

MPI runtimes and optional packages
---------------------------------

A serial-built PLUMED kernel can use the MD code's MPI through injection or a
runtime provider, without rebuilding PLUMED. Embedded thread-MPI is available
through the explicit `--tmpi-threads N` launcher. See the
[MPI wrapper documentation](src/mpiwrap/README.md) for runtime selection,
`plumed mkmpi`, compatibility requirements and current limitations.

This also allows one base `plumed` conda package with small optional
`plumed-mpi` variants for MPICH and Open MPI, rather than separate full PLUMED
builds. Our [conda recipes](conda/README.md) prototype and test that layout;
this does not change conda-forge's current packages.

Modules and plugins
-------------------

The source code is organized in modules, one for each directory of `src`.
The `module.type` file of each module says when it is compiled:

* `always`: the basic infrastructure (e.g. `core`, `tools`), always compiled;
* `default-on`: the common actions (e.g. `bias`, `colvar`, `function`,
  `generic`, `multicolvar`, `secondarystructure`), compiled unless they are
  disabled;
* `default-off`: specialized methods and modules that need external libraries
  (e.g. `opes`, `ves`, `funnel`, `sasa`, `pytorch`, `metatomic`), compiled only
  when they are requested.

Modules are selected when PLUMED is configured, with a `:`-separated list:

    ./configure --enable-modules=+opes:+ves   # add modules to the default ones
    ./configure --enable-modules=-isdb        # remove a default-on module
    ./configure --enable-modules=all          # all the modules
    ./configure --enable-modules=none:+opes   # only the always modules and opes

A module that needs an external library provides its actions only if PLUMED
is also configured with that library (see `./configure --help`, e.g.
`--enable-libtorch` for `pytorch`). The modules of an installed PLUMED are listed by
`plumed config show`, and `plumed config module opes` tells if a module is
enabled. The Python wheels (`pip install plumed`) and the CMake build for
Windows include only the `always` and `default-on` modules
(see [PLATFORMS.md](PLATFORMS.md)).

Plugins are compiled separately, against an installed PLUMED, and loaded at
run time with the `LOAD` action (`LOAD FILE=plugin.so` in the input, or
`plumed --load plugin.so` for the command line tools). This also works with the
Python wheels on Linux and macOS, but not on Windows. Plugins can be:

* your own actions: a C++ file with an action can be compiled with
  `plumed mklib MyCV.cpp`, or loaded directly with `LOAD FILE=MyCV.cpp`, which
  compiles it when PLUMED starts (a C++ compiler is needed in both cases);
* the plugins in the [plugins](plugins) directory, each with its own README and
  build instructions: [pycv](plugins/pycv/README.md) (collective variables and
  functions written in Python, installed with `pip`),
  [cudaCoord](plugins/cudaCoord/README.md) (coordination number on GPUs with
  CUDA) and [openaccPTM](plugins/openaccPTM/README.md) (collective variables on
  GPUs with OpenACC, including SAXS and SANS). cudaCoord and openaccPTM are
  compiled with `./configure` and `make` in their directory.

Python package and other packages
---------------------------------

The Python package, including a basic PLUMED kernel and the `plumed` command,
can be installed with

    python -m pip install plumed

Wheels are available for Linux (x86_64 and aarch64), macOS (Intel and Apple
Silicon) and Windows (x86_64). They are compiled without MPI libraries, OpenMP
and optional libraries; see [PLATFORMS.md](PLATFORMS.md) for the features
available on each platform and [python/README.rst](python/README.rst) for the
package itself.
PLUMED is also distributed by several package managers (see the badges above),
e.g. `conda install -c conda-forge plumed`.

Quick compilation instructions
------------------------------

Extensive installation instructions are in the [user documentation](https://www.plumed.org/doc).
Quick instructions:

    ./configure --prefix=$HOME/opt
    make
    make doc # optional
    make test # optional

User documentation can be found at `new-manual/site/index.html`.
Developer documentation can be found at `developer-doc/html/index.html`.
[Pre-compiled documentation](https://www.plumed.org/doc) is available online, so this is only required
if you are working with a modified version of the code!

In order to run PLUMED without installing it you should type `source sourceme.sh`. However,
we recommend installing PLUMED.
To install it in `$HOME/opt` (directory should be set during `./configure`):

    umask 022
    make install

Now you will be able to run plumed using e.g.

    plumed help

If you compiled your own documentation, paths to the installed documentation can be found with command `plumed info --user-doc`.

A sample modulefile with environment variable will be placed in
`$HOME/opt/lib/plumed/src/lib/modulefile`. This can be useful if you want to
install multiple PLUMED versions side by side and select them with env modules.

