#!/usr/bin/env bash
# Run against a provider-capable CLI and a native provider, without MPI_Init.
set -euo pipefail
if (($# != 3)); then
  echo 'usage: check-mpi-provider-selection.sh PROGRAM PROVIDER MPI_COMPILER' >&2
  exit 2
fi
program=$1
provider=$2
compiler=$3
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-provider-selection.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
mkdir "$TASK_DIR/paths with spaces" "$TASK_DIR/temporary builds"
export TMPDIR="$TASK_DIR/temporary builds"
cp "$provider" "$TASK_DIR/paths with spaces/provider.so"
PLUMED_MPI_WRAPPER="$TASK_DIR/paths with spaces/provider.so" \
  "$program" --no-mpi --has-real-mpi > "$TASK_DIR/library.log" 2>&1 || { cat "$TASK_DIR/library.log" >&2; exit 1; }

# Exercise the actual driver too: a Mach-O executable can be accepted by
# dlopen on macOS even though it is a compiler rather than a provider DSO.
PLUMED_MPI_WRAPPER="$compiler" \
  "$program" --no-mpi --has-real-mpi > "$TASK_DIR/direct-compiler.log" 2>&1 || { cat "$TASK_DIR/direct-compiler.log" >&2; exit 1; }

export MPI_TEST_COMPILER="$compiler"
export MPI_TEST_COMPILER_MARKER="$TASK_DIR/compiler-ran"
# This filename would expand if it were accidentally used as shell source.
selected_compiler="$TASK_DIR/paths with spaces/mpi\$(not-a-command)cc"
cat > "$selected_compiler" <<'COMPILER'
#!/usr/bin/env bash
printf 'compiler invoked\n' > "$MPI_TEST_COMPILER_MARKER"
exec "$MPI_TEST_COMPILER" "$@"
COMPILER
chmod +x "$selected_compiler"
PLUMED_MPI_WRAPPER="$selected_compiler" \
  "$program" --no-mpi --has-real-mpi > "$TASK_DIR/compiler.log" 2>&1 || { cat "$TASK_DIR/compiler.log" >&2; exit 1; }
if [[ ! -s "$MPI_TEST_COMPILER_MARKER" ]]; then
  echo 'provider selection did not invoke the selected compiler' >&2; exit 1
fi
printf '#!/usr/bin/env bash\nexit 7\n' > "$TASK_DIR/failing-compiler"
chmod +x "$TASK_DIR/failing-compiler"
# Retain the executable's upstream exception handling; require failure and a
# useful provider diagnostic without imposing a particular process exit code.
for selection in "$TASK_DIR/missing.so" relative/compiler "$TASK_DIR/failing-compiler"; do
  if PLUMED_MPI_WRAPPER="$selection" "$program" --no-mpi --has-real-mpi \
      > "$TASK_DIR/error.log" 2>&1; then
    echo "invalid provider selection succeeded: $selection" >&2; exit 1
  fi
  if ! grep -E 'MPI provider|PLUMED_MPI_WRAPPER' "$TASK_DIR/error.log" >/dev/null; then
    cat "$TASK_DIR/error.log" >&2; exit 1
  fi
done
cat > "$TASK_DIR/no-factory.cpp" <<'SOURCE'
#define PLUMED_MPIW_ABI_ONLY 1
#include "mpiwrap.h"
namespace PLMD {
namespace mpiwrap {
#ifdef TEST_BAD_FACTORY
extern "C" const plumed_mpiw_provider_descriptor *plumed_mpiw_provider_factory() { return 0; }
#else
extern "C" int unrelated_symbol() { return 0; }
#endif
}
}
SOURCE
case "$(uname -s)" in
  Darwin) shared_flags=(-dynamiclib) ;;
  *) shared_flags=(-shared) ;;
esac
"${CXX:-c++}" -std=c++98 -fPIC -iquote "$ROOT/src/mpiwrap" \
  "${shared_flags[@]}" "$TASK_DIR/no-factory.cpp" \
  -o "$TASK_DIR/no-factory.so"
chmod +x "$TASK_DIR/no-factory.so"
if PLUMED_MPI_WRAPPER="$TASK_DIR/no-factory.so" "$program" --no-mpi --has-real-mpi \
    > "$TASK_DIR/factory-error.log" 2>&1; then
  echo 'provider selection accepted a library without a factory' >&2; exit 1
fi
if ! grep -F 'no provider factory' "$TASK_DIR/factory-error.log" >/dev/null; then
  cat "$TASK_DIR/factory-error.log" >&2; exit 1
fi
"${CXX:-c++}" -std=c++98 -fPIC -DTEST_BAD_FACTORY -iquote "$ROOT/src/mpiwrap" \
  "${shared_flags[@]}" "$TASK_DIR/no-factory.cpp" -o "$TASK_DIR/bad-provider.so"
if PLUMED_MPI_WRAPPER="$TASK_DIR/bad-provider.so" "$program" --no-mpi --has-real-mpi \
    > "$TASK_DIR/abi-error.log" 2>&1; then
  echo 'provider selection accepted an invalid descriptor' >&2; exit 1
fi
if ! grep -F 'selected MPI provider is incompatible' "$TASK_DIR/abi-error.log" >/dev/null; then
  cat "$TASK_DIR/abi-error.log" >&2; exit 1
fi
shopt -s nullglob
remaining=("$TMPDIR"/plumed-mpi-provider.*)
if ((${#remaining[@]})); then
  echo 'provider selection left temporary build artifacts' >&2; exit 1
fi
echo 'environment provider library/compiler selection and cleanup passed'
