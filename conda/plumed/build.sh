#! /usr/bin/env bash

if [[ $(uname) == "Linux" ]]; then
  # STATIC_LIBS is a PLUMED specific option and is required on Linux for the following reason:
  # When using env modules the dependent libraries can be found through the
  # LD_LIBRARY_PATH or encoded configuring with -rpath.
  # Conda does not use LD_LIBRARY_PATH and it is thus necessary to suggest where libraries are.
  export STATIC_LIBS=-Wl,-rpath-link,$PREFIX/lib
  # The 2.28 sysroot keeps pthreads separate from libc. Library search is
  # disabled below, so both system libraries must be requested explicitly.
  export LIBS="-lpthread -lrt $LIBS"
fi

# we also store path so that software linking libplumedWrapper.a knows where libplumedKernel can be found.
export CPPFLAGS="-D__PLUMED_DEFAULT_KERNEL=$PREFIX/lib/libplumedKernel$SHLIB_EXT $CPPFLAGS"

make clean

# enable optimization
# and disable warnings, see below
export CXXFLAGS="${CXXFLAGS//-O2/-O3} -w"

if [[ $(uname) == "Darwin" ]]; then
  # see https://conda-forge.org/docs/maintainer/knowledge_base.html#newer-c-features-with-old-sdk
  CXXFLAGS="${CXXFLAGS} -D_LIBCPP_DISABLE_AVAILABILITY"

  # for some unexpected reason, when building on MacOS plumed 2.10 (likely because of C++17) # we need to pass the arguments below.
  # if we do not do so, linking plumed-runtime fails because the destructor of std::bad_function_call is not found
  export STATIC_LIBS="-Wl,-rpath,$PREFIX/lib -L$PREFIX/lib"
fi

# libraries are explicitly listed here due to --disable-libsearch
# (hdf5 for the HDF5 tables, files with extension .h5: the same library as h5py and other packages in the environment)
# Link the actual provider SONAME, not a switchable BLAS/LAPACK symlink.
if [[ $(uname) == Darwin ]]; then
  blas_library="$PREFIX/lib/libopenblas.0.dylib"
else
  blas_library="$PREFIX/lib/libopenblas.so.0"
fi
export LIBS="-lchemfiles -lfftw3 $blas_library -lz -lhdf5 $LIBS"

# python is disabled since it should be provided as a separate package
# --disable-libsearch forces to link only explicitely requested libraries
# --disable-static-patch avoid tests that are only required for static patches
# --disable-static-archive makes package smaller
# --disable-basic-warnings deactive the warning, compiling with the current compiler version (today is 28/10/2025) shows some false positives in ActionPilot
# --enable-libchemfiles reads more trajectory formats in the driver (--itraj)
./configure --prefix=$PREFIX --disable-mpi --disable-python --disable-libsearch --disable-static-patch --disable-static-archive --disable-basic-warnings --enable-libchemfiles || { cat config.log; exit 1; }

make -j${CPU_COUNT}
make install
# Runtime compilation is supplied by the optional plumed-mklib package.
rm -f "$PREFIX/lib/plumed/scripts/mklib.sh" "$PREFIX/lib/plumed/plumed-mklib"
