<!--
description: Using arbitrary machine learning models as collective variables
authors: Guillaume Fraux
reference:
-->

# Overview

This module implements the interface between PLUMED and [metatomic], allowing to
use arbitrary machine learning models as collective variables. These machine
learning models are defined using custom Python code — following the [metatomic]
interface — and then either exported to [TorchScript], or loaded directly from
Python scripts. The model is then loaded inside PLUMED and executed during the
simulation.

PLUMED uses the C++ API of metatomic-core to load and execute models, through
metatomic plugins: the plugin of `metatomic-torch` runs TorchScript models
without starting Python, and the Python plugin of metatomic-core runs models
defined in Python scripts. See the documentation of the [METATOMIC](METATOMIC.md)
action for how to select plugins.

# Installation

This module requires the C++ metatomic (from metatomic-core) and metatensor
libraries. TorchScript models also need the plugin of `metatomic-torch`.

## Getting the dependencies

The easiest way to get all dependencies is to install the Python packages with
`pip`, and use the libraries they contain. PLUMED then finds the Python plugin
of the `metatomic` package by itself. Give the path of the `metatomic-torch`
plugin (`python -c "import metatomic_torch; print(metatomic_torch.utils.plugin_path)"`)
to the `PLUGIN_LIBRARY` keyword of [METATOMIC](METATOMIC.md) to run TorchScript
models.

metatomic-core is not yet released on PyPI, and needs to be installed from the
sources of [metatomic], together with `metatomic-torch` for TorchScript models.
The script `.ci/install.metatomic` shows how this is done for PLUMED's tests.
Building from source requires Python 3.11 or newer, a Rust compiler (version
1.88 or newer) and CMake.

```bash
# you can then get the compiler and linker flags using the script at
# src/metatomic/flags-from-python.py
CPPFLAGS=$(python src/metatomic/flags-from-python.py --cppflags)
LDFLAGS=$(python src/metatomic/flags-from-python.py --ldflags)
```

## Building plumed with metatomic

Once you installed all dependencies, you can now configure PLUMED:

```bash
./configure --enable-libmetatomic \
            --enable-modules=+metatomic \
            LDFLAGS="$LDFLAGS" CPPFLAGS="$CPPFLAGS"
```

!!! note

    Pay close attention to the output, it should contain a line about
    `checking libmetatomic` ending with `...yes`. If this is not the case,
    you'll get a warning about `cannot enable __PLUMED_HAS_LIBMETATOMIC`. You
    should then check `config.log` to know more about what's going on and why
    the library can't be found.

## Running models

Models run on CPU only (`DEVICE=cpu`).

To bias a model output, the model must compute its gradients with respect to
`positions` and `strain` explicitly, since metatomic only computes gradients of
energies with automatic differentiation. Outputs without these gradients can be
used and printed, but applying a bias to them stops PLUMED with an error.

Models defined in Python scripts run in the Python interpreter started by the
Python plugin, which is the one given by the `METATOMIC_PYTHON` environment
variable, or the active virtual or conda environment, or the first `python` in
`PATH`. If the same PLUMED run also uses [pycv](https://github.com/plumed/plumed2/tree/master/plugins/pycv),
build pycv with the same Python installation, since a process can only contain
one Python interpreter.


[TorchScript]: https://pytorch.org/docs/stable/jit.html
[metatomic]: https://docs.metatensor.org/metatomic/
