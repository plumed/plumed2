#!/bin/bash

env | sort

# GB: install xdrfile library
wget http://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
tar xzf xdrfile-1.1.4.tar.gz
cd xdrfile-1.1.4
./configure --prefix=$PREFIX --enable-shared
make
make install
cd ../

# TG: The "disabled" features are workaround for possible
#     conda+configure bugs in library search: building is ok but
#     linking with the .so doesn't find them (in
#     conda-forge). Possibly the LD path needs tweaks.

# TODO: re-enable them and see. Also to do: install docs?

./configure --prefix=$PREFIX --enable-shared --disable-python --disable-zlib --disable-external-lapack --disable-external-blas LDFLAGS=-L$PREFIX/lib
make -j4
make install

cd python
make pip
export plumed_default_kernel=$PREFIX/lib/libplumedKernel$SHLIB_EXT
$PYTHON -m pip install .

