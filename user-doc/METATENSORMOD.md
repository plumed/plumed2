\page METATENSORMOD Metatensor

<!--
description: Using arbitrary machine learning models as collective variables
authors: Guillaume Fraux
reference:
-->

# Overview

This module implements the interface between PLUMED and [metatensor], allowing
to use arbitrary machine learning models as collective variables. These machine
learning models are defined using custom Python code — following the [metatensor
atomistic models][mts_models] interface — and then exported to [TorchScript].
The exported model is then loaded inside PLUMED and executed during the
simulation.


# Installation

This module requires two main dependencies: the C++ torch library (i.e.
`libtorch`); and the C++ metatensor_torch library. There are multiple ways of
installing both libraries, which are discussed below.

## Installing the libraries through Python's package manager (`pip`)

The easiest way to get all dependencies on your system is to download the
pre-built Python wheels with `pip`. This is the same set of wheels you will need
to define custom models.

```bash
pip install "metatensor-torch ==0.5.5"  # change this version to get newer releases

# optional: get the other metatensor tools to define models (these are only usable from Python).
pip install metatensor-operations metatensor-learn

# export the location to all the libraries:
TORCH_CMAKE_PREFIX=$(python -c "import torch; print(torch.utils.cmake_prefix_path)")
TORCH_PREFIX=$(cd "$TORCH_CMAKE_PREFIX/../.." && pwd)

METATENSOR_CMAKE_PREFIX=$(python -c "import metatensor; print(metatensor.utils.cmake_prefix_path)")
METATENSOR_PREFIX=$(cd "$METATENSOR_CMAKE_PREFIX/../.." && pwd)

METATENSOR_TORCH_CMAKE_PREFIX=$(python -c "import metatensor.torch; print(metatensor.torch.utils.cmake_prefix_path)")
METATENSOR_TORCH_PREFIX=$(cd "$METATENSOR_TORCH_CMAKE_PREFIX/../.." && pwd)

# The torch library installed by pip uses a pre-cxx11 ABI
TORCH_CPPFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"
```

That's it, you can now jump to [the last part](#building-plumed-with-metatensor)
of the installation instructions.

## Using pre-built libraries

If you only want to use existing models, you can download pre-built versions of
the libraries and build PLUMED against these. First, you'll need to download
libtorch (see also \ref installation-libtorch for other instructions on
installing a pre-built libtorch):

```bash
# Download torch 2.2.2 for x86_64 (Intel) Linux.
#
# Variations of this link for other operating systems (macOS, Windows), CPU
# architecture (Apple Silicon, arm64), CUDA versions, and newer versions of
# libtorch can be found at https://pytorch.org/get-started/locally/

wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.2.2%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-2.2.2+cpu.zip

# alternatively if you have a CUDA-enabled GPU, you can use the corresponding
# pre-built library (here for CUDA 12.1):
wget https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-2.2.2%2Bcu121.zip
unzip libtorch-cxx11-abi-shared-with-deps-2.2.2+cu121.zip

# Make the location of libtorch visible
TORCH_PREFIX=$(pwd)/libtorch

# if you are using a library with pre-cxx11 ABI, you need an extra flag:
TORCH_CPPFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"
```

Once you acquire libtorch, it is time to build metatensor and metatensor_torch
from sources. There is currently no standalone pre-built library for these
(although you can use the pre-built version that comes with `pip`). For this,
you'll need a rust compiler on your system, which you can get with
[rustup](https://rustup.rs/) or any other method at your convenience.

```bash
# patch a bug from torch's MKL detection in CMake
cd <PLUMED/DIR>
./src/metatensor/patch-torch.sh "$TORCH_PREFIX"

cd <SOME/PLACE/WHERE/TO/PUT/METATENSOR/SOURCES>

# define a location where metatensor should be installed
METATENSOR_PREFIX=<...>

METATENSOR_TORCH_PREFIX="$METATENSOR_PREFIX"

git clone https://github.com/lab-cosmo/metatensor
# or a more recent release of metatensor-torch
git checkout metatensor-torch-v0.5.5
cd metatensor

mkdir build && cd build
cmake -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX="$METATENSOR_PREFIX" \
      -DCMAKE_PREFIX_PATH="$TORCH_PREFIX" \
      -DBUILD_METATENSOR_TORCH=ON \
      -DMETATENSOR_INSTALL_BOTH_STATIC_SHARED=OFF \
      ..

cmake --build . --target install --parallel
```

## Building plumed with metatensor

Once you installed all dependencies with one of the methods above, you can now
configure PLUMED:

```bash
# set include search path for the compilers
TORCH_INCLUDES="-I$TORCH_PREFIX/include -I$TORCH_PREFIX/include/torch/csrc/api/include"
CPPFLAGS="$TORCH_INCLUDES $TORCH_CPPFLAGS -I$METATENSOR_PREFIX/include -I$METATENSOR_TORCH_PREFIX/include $CPPFLAGS"

# set library search path for the linker
LDFLAGS="-L$TORCH_PREFIX/lib -L$METATENSOR_PREFIX/lib -L$METATENSOR_TORCH_PREFIX/lib $LDFLAGS"

# set the rpath to make sure plumed executable will be able to find the right libraries
LDFLAGS="$LDFLAGS -Wl,-rpath,$TORCH_PREFIX/lib"
LDFLAGS="$LDFLAGS -Wl,-rpath,$METATENSOR_PREFIX/lib -Wl,-rpath,$METATENSOR_TORCH_PREFIX/lib"

# If you are running on Linux, force the use of rpath instead of runpath
# (we rely on the rpath to find dependencies of dependencies)
LDFLAGS="$LDFLAGS -Wl,--disable-new-dtags"

# configure PLUMED
./configure --enable-libtorch --enable-metatensor --enable-modules=+metatensor \
    LDFLAGS="$LDFLAGS" CPPFLAGS="$CPPFLAGS"
```

Pay close attention to the output, it should contain **both** a line about
`checking libtorch` and a line about `checking metatensor`, both ending with
`...yes`. If this is not the case, you'll get a warning about `cannot enable
__PLUMED_HAS_LIBTORCH` or `cannot enable __PLUMED_HAS_METATENSOR`. If you get
any of these warnings, you should check `config.log` to know more about what's
going on and why these libraries can't be found.

# Module Contents

This module defines the following actions:

@METATENSORMOD_COLVAR@




[TorchScript]: https://pytorch.org/docs/stable/jit.html
[metatensor]: https://docs.metatensor.org/latest/index.html
[mts_models]: https://docs.metatensor.org/latest/atomistic/index.html
