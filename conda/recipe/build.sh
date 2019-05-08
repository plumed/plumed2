#!/bin/bash

env | sort

if [[ $target_platform == linux* ]]; then
    condaldflags="$BUILD_PREFIX/lib/libz.a $BUILD_PREFIX/lib/libz.so.1 $BUILD_PREFIX/lib/libxdrfile.a $BUILD_PREFIX/lib/libxdrfile.so"
elif [[ $target_platform == osx* ]]; then
    # Possibly not needed.
    condaldflags="$BUILD_PREFIX/lib/libz.a $BUILD_PREFIX/lib/libz.dylib $BUILD_PREFIX/lib/libxdrfile.a $BUILD_PREFIX/lib/libxdrfile.dylib"
else
    echo "Unexpected platform $target_platform"
fi



# GB: install xdrfile library
wget http://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
tar xzf xdrfile-1.1.4.tar.gz
cd xdrfile-1.1.4
./configure --prefix=$PREFIX --enable-shared
make
make install
cd ../



# Thanks to https://github.com/intbio/plumed-conda/blob/master/plumed2_v2.5.0/build.sh

./configure --prefix=$PREFIX --enable-shared --disable-python --disable-external-lapack --disable-external-blas LDFLAGS="$condaldflags"
make -j4
make install

cd python
make pip
export plumed_default_kernel=$PREFIX/lib/libplumedKernel$SHLIB_EXT
$PYTHON -m pip install .

