#!/bin/bash

env | sort

# GB: install xdrfile library
if true; then
    wget http://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
    tar xzf xdrfile-1.1.4.tar.gz
    cd xdrfile-1.1.4
    ./configure --prefix=$PREFIX --enable-shared
    make
    make install
    cd ../
fi

# TODO: install docs?

# python wrapper is installed with pip
# we temporarily use internal lapack/blas (should probably be fixed)

if test -n "$MACOSX_DEPLOYMENT_TARGET" ; then
  opt=""
else
# STATIC_LIBS is required on Linux for the following reason:
# When using env modules the dependent libraries can be found through the
# LD_LIBRARY_PATH or encoded configuring with -rpath.
# Conda does not use LD_LIBRARY_PATH and it is thus necessary to suggest where libraries are.
  opt=STATIC_LIBS=-Wl,-rpath-link,$PREFIX/lib
# Needed due to asmjit. See
# See https://github.com/Bumblebee-Project/Bumblebee/issues/76
  export LIBS="$LIBS -lrt"
fi

./configure --prefix=$PREFIX --enable-shared --disable-python --disable-external-lapack --disable-external-blas --enable-asmjit $opt

make -j4
make install

