#!/usr/bin/env bash
set -euo pipefail
export PLUMED_CONDA_BUILD=1
# The generic extension route (extensions/build.py), as for the other modules:
# pycv (toolchain python) is compiled for the Python of this variant, with the
# pybind11 of the host environment, linked to its libpython and registered with
# the variant python<version>.
./configure --with-extension-base="$PREFIX" --enable-extension-modules=pycv \
  --prefix="$PREFIX" PYTHON_BIN="$PYTHON"
make
make install
