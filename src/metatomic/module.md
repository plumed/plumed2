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

PLUMED uses the C++ API of metatomic-core to load and execute models.
TorchScript models exported with `metatomic-torch` run through the Python plugin
of metatomic-core, which starts the Python interpreter of the active environment
inside PLUMED. Models can also be loaded through other metatomic plugins, see
the documentation of the [METATOMIC](METATOMIC.md) action.


# Installation

This module requires the C++ metatomic (from metatomic-core) and metatensor
libraries. To use TorchScript models, the Python environment found by the
metatomic Python plugin also needs `metatomic-torch`.

## Getting the dependencies

The easiest way to get all dependencies is to install the Python packages with
`pip`, and use the libraries they contain. This ensures that PLUMED and the
Python plugin share the same copy of `libmetatomic`.

metatomic-core is not yet released on PyPI, and needs to be installed from the
sources of [metatomic], together with `metatomic-torch` for TorchScript models.
The script `.ci/install.metatomic` shows how this is done for PLUMED's tests,
with the metatomic commit they use. Building from source requires Python 3.11
or newer, a Rust compiler (version 1.88 or newer) and CMake.

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

At runtime, PLUMED registers the Python plugin installed next to
`libmetatomic`. The plugin uses the Python interpreter given by the
`METATOMIC_PYTHON` environment variable, or the active virtual or conda
environment, or the first `python` in `PATH`. This interpreter must be able to
import `metatomic` (and `metatomic_torch` for TorchScript models), from the same
installation PLUMED was built with.

This PLUMED backend supports `DEVICE=cpu` only. This is a restriction of the
backend, not of `LegacyModel`.

To bias an output, the model must provide both position and strain gradients.
At the pinned core version, the standard `feature` output has no autograd
support, for either per-system or per-atom samples. The model must declare and
return both gradients explicitly. Features without these gradients can be
printed, but applying a bias produces an error.

Core consistency checks currently reject gradients on `feature`. Leave
`CHECK_CONSISTENCY` off when biasing these outputs. The per-atom regression test
checks forces and all nine strain derivatives numerically; this does not replace
a future core fix to accept feature gradients.

If you also use pycv, build and run both plugins with the same Python
installation. Set `METATOMIC_PYTHON` to that interpreter before starting PLUMED:

```bash
export METATOMIC_PYTHON="$(python -c 'import sys; print(sys.executable)')"
# In plugins/pycv, use this same interpreter to build and test:
CMAKE_ARGS="-DPython_EXECUTABLE=$METATOMIC_PYTHON" make check PYTHON="$METATOMIC_PYTHON"
```

If pycv initializes Python inside a virtual environment, also make that
environment's packages visible through `PYTHONPATH`. `METATOMIC_PYTHON` selects
the interpreter only when the metatomic plugin initializes Python itself.

The CI workflow uses one Python 3.11 installation for dependency builds, model
exports, pycv, and the metatomic Python plugin.


[TorchScript]: https://pytorch.org/docs/stable/jit.html
[metatomic]: https://docs.metatensor.org/metatomic/
