#!/bin/bash

if [[ $(uname) == "Linux" ]]; then
# STATIC_LIBS is a PLUMED specific option and is required on Linux for the following reason:
# When using env modules the dependent libraries can be found through the
# LD_LIBRARY_PATH or encoded configuring with -rpath.
# Conda does not use LD_LIBRARY_PATH and it is thus necessary to suggest where libraries are.
  export STATIC_LIBS=-Wl,-rpath-link,$PREFIX/lib
# -lrt is required to link clock_gettime
  export LIBS="-lrt $LIBS"
fi

# we also store path so that software linking libplumedWrapper.a knows where libplumedKernel can be found.
export CPPFLAGS="-D__PLUMED_DEFAULT_KERNEL=$PREFIX/lib/libplumedKernel$SHLIB_EXT $CPPFLAGS"

# enable optimization
export CXXFLAGS="${CXXFLAGS//-O2/-O3}"

# libraries are explicitly listed here due to --disable-libsearch
export LIBS="-lfftw3 -lgsl -lgslcblas -llapack -lblas -lxdrfile -lz $LIBS"

# python is disabled since it should be provided as a separate package
# --disable-libsearch forces to link only explicitely requested libraries
# --disable-static-patch avoid tests that are only required for static patches
# --disable-static-archive makes package smaller
# --enable-asmjit enables bundled asmjit implementation
./configure --prefix=$PREFIX --disable-python --disable-libsearch --disable-static-patch --disable-static-archive --enable-asmjit

make -j${CPU_COUNT}
make install

