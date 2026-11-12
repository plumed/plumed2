#! /usr/bin/env bash

cd python
# Local source trees can contain a wrapper compiled for a different prefix.
rm -rf build
make pip
export PLUMED_BUILD_KERNEL=0
export plumed_default_kernel=$PREFIX/lib/libplumedKernel$SHLIB_EXT
#export plumed_disable_rtld_deepbind=yes

$PYTHON -m pip install . --no-deps -vv

