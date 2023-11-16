#! /usr/bin/env bash

cd python
make pip
export plumed_default_kernel=$PREFIX/lib/libplumedKernel$SHLIB_EXT
export plumed_disable_rtld_deepbind=yes

$PYTHON -m pip install . --no-deps -vv

