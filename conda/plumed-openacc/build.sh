#!/usr/bin/env bash
set -euo pipefail
export PLUMED_AUTOLOAD=0
export PLUMED_KERNEL="$PREFIX/lib/libplumedKernel.so"
: "${NVHPC_HOME:?set NVHPC_HOME to the isolated NVIDIA HPC SDK installation}"
sdk="$NVHPC_HOME/Linux_$(uname -m)/$nvhpc_version"
nvcxx="$sdk/compilers/bin/nvc++"
"$nvcxx" -V
mkdir -p "$SRC_DIR/nvhpc-licenses"
cp "$sdk/compilers/license/"* "$SRC_DIR/nvhpc-licenses/"
cp "$sdk/cuda/$openacc_cuda_version/EULA.txt" "$SRC_DIR/nvhpc-licenses/CUDA-EULA.txt"
python extensions/register.py --plumed "$PREFIX/bin/plumed" --module openacc_ptm \
  --library ../../plumedOpenACC.so --variant "nvhpc${nvhpc_version}_cuda${openacc_cuda_version}" \
  --features OPENACC EXTERNAL_LAPACK --output "$SRC_DIR/openacc-registration"
cp "$SRC_DIR/openacc-registration/module_identity.cpp" "$SRC_DIR/plugins/openaccPTM/src/ModuleIdentity.cpp"
# NVHPC needs explicit GNU compiler locations: conda activation exports names,
# not an installation root containing unprefixed gcc/g++ binaries. Generate a
# private profile; never rewrite the externally supplied SDK configuration.
"$sdk/compilers/bin/makelocalrc" -d "$SRC_DIR/openacc-registration" \
  -x "$sdk/compilers/bin" -gcc "$(command -v "$CC")" \
  -gpp "$(command -v "$CXX")" -g77 "$(command -v gfortran)"
export NVLOCALRC="$SRC_DIR/openacc-registration/localrc"
test -s "$NVLOCALRC"
sysroot=$("$CC" -print-sysroot)
test -d "$sysroot/usr/include"
# NVHPC otherwise searches host glibc libraries even with conda's GCC profile.
# Put the conda sysroot first and apply it to GNU ld's library scripts as well.
sysroot_link="-L$sysroot/usr/lib64 -L$sysroot/lib64 -L$sysroot/usr/lib -L$sysroot/lib -Xlinker --sysroot=$sysroot"
policy="$SRC_DIR/extensions/conda_config.py"
cd plugins/openaccPTM
python "$policy" --module openacc_ptm -- bash configure
make clean
# Include real GPU machine code for all supported architectures, and a serial
# host target so ordinary CI can test the same installed artifact without a GPU.
# Keep the vendor runtime private: the base uses GCC's OpenMP runtime.
# Preserve the base CPPFLAGS, including MPIW and other public-header ABI macros.
python "$policy" --module openacc_ptm -- make -j"${CPU_COUNT:-4}" NVCXX="$nvcxx" ACCTARGET=gpu,host \
  CXXFLAGS="-O2 -std=c++17 -fPIC -tp=px -gpu=ccall,cuda${openacc_cuda_version}" \
  NVLDSHARED="$nvcxx -shared $sysroot_link -tp=px -gpu=ccall,cuda${openacc_cuda_version} -Xlinker --no-undefined" \
  DYNAMIC_LIBS="-L$PREFIX/lib $PREFIX/lib/libopenblas.so.0"
readelf -d plumedOpenACC.so
python - <<'PYCODE'
import re
import subprocess
versions = subprocess.check_output(["readelf", "--version-info", "plumedOpenACC.so"], text=True)
required = max((int(a), int(b)) for a, b in re.findall(r"GLIBC_(\d+)\.(\d+)", versions))
print("OpenACC minimum glibc:", ".".join(map(str, required)), flush=True)
assert required <= (2, 28), "OpenACC linked against host glibc instead of conda sysroot"
PYCODE
# Do not ship a library that depends on an SDK installation on the build host.
if readelf -d plumedOpenACC.so | grep NEEDED | grep -E 'lib(acc|nv|pg)'; then
  echo 'Unexpected shared NVIDIA runtime dependency' >&2
  exit 1
fi
install -m 755 plumedOpenACC.so "$PREFIX/lib/"
mkdir -p "$PREFIX/lib/plumed/plugins"
cp "$SRC_DIR/openacc-registration/openacc_ptm.plumed" "$PREFIX/lib/plumed/plugins/"

python "$policy" --module openacc_ptm --check-features --prefix "$PREFIX"
