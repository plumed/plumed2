#!/usr/bin/env bash
# Exercise the real mkmpi discovery path with a deterministic mock compiler.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiw-discovery.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/scripts" "$tmp/src/config" "$tmp/src/wrapper"
cp "$root/scripts/mkmpi.sh" "$tmp/scripts/"
cp "$root/src/wrapper/Plumed.h" "$tmp/src/wrapper/"
printf 'mpi_provider_link_flags="-shared"\n' > "$tmp/src/config/compile_options.sh"
case $(uname -s) in
  Darwin) shared=-dynamiclib; extension=dylib ;;
  *) shared=-shared; extension=so ;;
esac
"${CXX:-c++}" -std=c++17 -fPIC "$shared" -DMPIW_TEST_PROVIDER \
  -I"$root/src" -iquote "$root/src/mpiwrap" \
  -x c++ "$root/regtest/basic/rt-mpiwrap-api/test-provider-discovery.cpp.in" \
  -o "$tmp/provider.so"
"${CXX:-c++}" -std=c++17 -pthread -D__PLUMED_MPIWRAP -D__PLUMED_HAS_DLOPEN \
  -I"$root/src" -iquote "$root/src/mpiwrap" \
  -x c++ "$root/src/mpiwrap/mpiwrap.cpp" \
  -x c++ "$root/src/mpiwrap/loader.cpp" \
  -x c++ "$root/src/mpiwrap/fallbacks.cpp" \
  -x c++ "$root/src/mpiwrap/launcher.cpp" \
  -x c++ "$root/src/mpiwrap/tmpi.cpp" \
  -x c++ "$root/regtest/basic/rt-mpiwrap-api/test-provider-discovery.cpp.in" -ldl -o "$tmp/test"
cat > "$tmp/mock-mpicc" <<'COMPILER'
#!/usr/bin/env bash
set -euo pipefail
printf 'compile\n' >> "$MPIW_TEST_ROOT/compiler.log"
while (($#)); do
  if [[ $1 == -o ]]; then cp "$MPIW_TEST_ROOT/provider.so" "$2"; exit; fi
  shift
done
exit 1
COMPILER
chmod +x "$tmp/mock-mpicc"
unset PLUMED_MPI_WRAPPER PLUMED_MPI_WRAPPER_FORCE
unset MPIW_TEST_KERNEL_PATH MPIW_TEST_PROGRAM_NAME MPIW_TEST_REJECT_DESCRIPTOR
export MPIW_TEST_ROOT="$tmp" MPICC="$tmp/mock-mpicc"
"$tmp/test" success
[[ $(wc -l < "$tmp/compiler.log") -eq 1 ]]
# Explicitly failing discovery must leave the registry available for later use.
export MPICC=/usr/bin/false
"$tmp/test" failure
[[ $(wc -l < "$tmp/compiler.log") -eq 2 ]]
export MPICC="$tmp/mock-mpicc"
# A relocated kernel finds its sibling provider, following the kernel symlink.
# Ordinary constants/inquiries in the fixture must leave this provider unselected.
mkdir -p "$tmp/relocated lib" "$tmp/kernel links"
touch "$tmp/relocated lib/libplumedKernel.$extension"
ln -s "../relocated lib/libplumedKernel.$extension" "$tmp/kernel links/kernel.$extension"
export MPIW_TEST_KERNEL_PATH="$tmp/kernel links/kernel.$extension"
cp "$tmp/provider.so" "$tmp/relocated lib/libplumed-mpi-provider.$extension"
"$tmp/test" installed
[[ $(wc -l < "$tmp/compiler.log") -eq 2 ]]
MPIW_TEST_REJECT_DESCRIPTOR=1 "$tmp/test" installed-invalid
[[ $(wc -l < "$tmp/compiler.log") -eq 2 ]]
# Custom program names use the same provider suffix without rebuilding the core.
mv "$tmp/relocated lib/libplumed-mpi-provider.$extension" \
  "$tmp/relocated lib/libplumed_custom-mpi-provider.$extension"
MPIW_TEST_PROGRAM_NAME=plumed_custom "$tmp/test" installed
[[ $(wc -l < "$tmp/compiler.log") -eq 2 ]]
# Missing siblings retain the existing compilation fallback.
rm "$tmp/relocated lib/libplumed_custom-mpi-provider.$extension"
"$tmp/test" success
[[ $(wc -l < "$tmp/compiler.log") -eq 3 ]]
# A present but unusable provider is an error, never a compiler fallback.
printf 'invalid library\n' > "$tmp/relocated lib/libplumed-mpi-provider.$extension"
"$tmp/test" installed-invalid
[[ $(wc -l < "$tmp/compiler.log") -eq 3 ]]
# Injection and explicit environment selection bypass even an invalid sibling.
"$tmp/test" injection
"$tmp/test" environment
[[ $(wc -l < "$tmp/compiler.log") -eq 3 ]]
rm "$tmp/relocated lib/libplumed-mpi-provider.$extension"
ln -s missing-provider "$tmp/relocated lib/libplumed-mpi-provider.$extension"
"$tmp/test" installed-invalid
[[ $(wc -l < "$tmp/compiler.log") -eq 3 ]]
echo 'native communicator provider discovery and failure handling passed'
