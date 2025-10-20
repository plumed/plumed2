<!--
description: Using arbitrary machine learning models as collective variables
authors: Guillaume Fraux
reference:
-->

# Overview

This module implements the interface between PLUMED and [metatomic], allowing to
use arbitrary machine learning models as collective variables. These machine
learning models are defined using custom Python code — following the [metatomic]
interface — and then exported to [TorchScript]. The exported model is then
loaded inside PLUMED and executed during the simulation.


# Installation

This module requires two main dependencies: the C++ torch library (i.e.
`libtorch`); and the C++ metatomic_torch library. There are multiple ways of
installing both libraries, which are discussed below.

## Getting the dependencies

### Using Python's package manager (`pip`)

The easiest way to get all dependencies on your system is to download the
pre-built Python wheels with `pip`. This is the same set of wheels you will need
to define custom models.

```bash
# change this version to get newer releases
pip install "metatomic-torch ==0.1.2"

# you can then get the compiler and linker flags using the script at
# src/metatomic/flags-from-python.py
CPPFLAGS=$(python src/metatomic/flags-from-python.py --cppflags)
LDFLAGS=$(python src/metatomic/flags-from-python.py --ldflags)
```

That's it, you can now jump to [the last part](#building-plumed-with-metatensor)
of the installation instructions.

### Using conda

You can also use conda with the [conda-forge] channel to install all
dependencies:

```bash
conda install -c conda-forge libmetatomic-torch

# setup compiler flags:
CPPFLAGS="-I$CONDA_PREFIX/include $CPPFLAGS"
CPPFLAGS="-I$CONDA_PREFIX/include/torch/csrc/api/include $CPPFLAGS"

LDFLAGS="-L$CONDA_PREFIX/lib $LDFLAGS"
LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib $LDFLAGS"

# if you are using linux, force the use of rpath instead of runpath
# (we rely on the rpath to find dependencies of dependencies)
LDFLAGS="-Wl,--disable-new-dtags"
```

That's it, you can now jump to [the last part](#building-plumed-with-metatensor)
of the installation instructions.

### From sources

You can also compile all required libraries from source. To this end, please
follow the corresponding instructions:

- **torch**: <https://pytorch.org/get-started/locally/>, for the C++ language
- **metatensor**: <https://docs.metatensor.org/latest/installation.html#install-c>
- **metatensor-torch**: <https://docs.metatensor.org/latest/installation.html#install-torch-cxx>
- **metatomic-torch**: <https://docs.metatensor.org/metatomic/latest/installation.html#install-torch-cxx>

Once everything is compiled and installed, you should set the corresponding
preprocessor flags in `CPPFLAGS` and the corresponding linker flags in
`LDFLAGS`.

## Building plumed with metatensor

Once you installed all dependencies with one of the methods above, you can now
configure PLUMED:

```bash
./configure --enable-libtorch \
            --enable-libmetatomic \
            --enable-modules=+metatomic \
            LDFLAGS="$LDFLAGS" CPPFLAGS="$CPPFLAGS"
```

!!! note

    Pay close attention to the output, it should contain **both** a line about
    `checking libtorch` and a line about `checking libmetatomic`, both ending with
    `...yes`. If this is not the case, you'll get a warning about `cannot enable
    __PLUMED_HAS_LIBTORCH` or `cannot enable __PLUMED_HAS_LIBMETATOMIC`.

    If you get any of these warnings, you should check `config.log` to know more
    about what's going on and why these libraries can't be found.


[TorchScript]: https://pytorch.org/docs/stable/jit.html
[metatomic]: https://docs.metatensor.org/metatomic/
[conda-forge]: https://conda-forge.org/
