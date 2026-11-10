#!/usr/bin/env bash
set -euo pipefail
export PLUMED_AUTOLOAD=0
export PLUMED_KERNEL="$PREFIX/lib/libplumedKernel.so"
: "${NVHPC_HOME:?set NVHPC_HOME to the isolated NVIDIA HPC SDK installation}"
sdk="$NVHPC_HOME/Linux_$(uname -m)/$nvhpc_version"
nvcxx="$sdk/compilers/bin/nvc++"
"$nvcxx" -V
mkdir -p "$SRC_DIR/nvhpc-licenses"
cp "$sdk/compilers/license/"LICENSE*.txt "$SRC_DIR/nvhpc-licenses/"
python extensions/register.py --plumed "$PREFIX/bin/plumed" --module openacc_ptm \
  --library ../../plumedOpenACC.so --variant "nvhpc${nvhpc_version}_cuda${openacc_cuda_version}" \
  --features OPENACC --output "$SRC_DIR/openacc-registration"
"$CXX" $CXXFLAGS $CPPFLAGS -std=c++17 -fPIC -D__PLUMED_HAS_OPENACC=1 \
  -I"$PREFIX/include/plumed" -c "$SRC_DIR/openacc-registration/module_identity.cpp" \
  -o "$SRC_DIR/openacc-registration/module_identity.o"
cd plugins/openaccPTM
bash configure
make clean
# Include real GPU machine code for all supported architectures, and a serial
# host target so ordinary CI can test the same installed artifact without a GPU.
# Keep the vendor runtime private: the base uses GCC's OpenMP runtime.
make -j"${CPU_COUNT:-4}" NVCXX="$nvcxx" ACCTARGET=gpu,host \
  CPPFLAGS="-I$PREFIX/include" \
  CXXFLAGS="-O2 -std=c++17 -fPIC --gcc-toolchain=$CC -gpu=ccall,cuda${openacc_cuda_version}" \
  NVLDSHARED="$nvcxx -shared --gcc-toolchain=$CC -gpu=ccall,cuda${openacc_cuda_version} -Xlinker --no-undefined" \
  EXTRA_OBJS="$SRC_DIR/openacc-registration/module_identity.o" \
  DYNAMIC_LIBS="-L$PREFIX/lib"
readelf -d plumedOpenACC.so
# Do not ship a library that depends on an SDK installation on the build host.
if readelf -d plumedOpenACC.so | grep NEEDED | grep -E 'lib(acc|nv|pg)'; then
  echo 'Unexpected shared NVIDIA runtime dependency' >&2
  exit 1
fi
install -m 755 plumedOpenACC.so "$PREFIX/lib/"
mkdir -p "$PREFIX/lib/plumed/plugins"
cp "$SRC_DIR/openacc-registration/openacc_ptm.plumed" "$PREFIX/lib/plumed/plugins/"
