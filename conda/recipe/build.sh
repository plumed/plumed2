#!/bin/bash

env | sort

# TG: The "disabled" features are workaround for possible conda+configure bugs in library
#     search: building is ok but linking with the .so doesn't find them.

./configure --prefix=$PREFIX --enable-shared --disable-python --disable-zlib --disable-external-lapack --disable-external-blas
make -j4
make install

cd python
make pip
export plumed_default_kernel=$PREFIX/lib
$PYTHON -m pip install 
