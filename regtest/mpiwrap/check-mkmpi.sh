#!/usr/bin/env bash
# Exercise argument boundaries, installation independence, and failed publication.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-mkmpi-test.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
export PLUMED_ROOT="$TASK_DIR/installed tree"
mkdir -p "$PLUMED_ROOT/src/config" "$PLUMED_ROOT/include/plumed/wrapper" "$TASK_DIR/compiler dir" "$TASK_DIR/output dir"
printf '/* public header */\n' > "$PLUMED_ROOT/include/plumed/wrapper/Plumed.h"
export PLUMED_INCLUDEDIR="$PLUMED_ROOT/include"
export PLUMED_PROGRAM_NAME=plumed
printf "mpi_provider_compile_flags='-fPIC'\nmpi_provider_link_flags='-shared'\n" > "$PLUMED_ROOT/src/config/compile_options.sh"
export ARGUMENT_LOG="$TASK_DIR/arguments"
cat > "$TASK_DIR/compiler dir/mpicc" <<'COMPILER'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$@" > "$ARGUMENT_LOG"
output=
while (($#)); do
  if [[ "$1" == -o ]]; then output=$2; shift; fi
  shift
done
if [[ "${FAIL_COMPILE:-no}" == yes ]]; then exit 2; fi
printf 'new provider\n' > "$output"
COMPILER
chmod +x "$TASK_DIR/compiler dir/mpicc"
output="$TASK_DIR/output dir/provider.so"
printf 'previous provider\n' > "$output"
if FAIL_COMPILE=yes "$ROOT/scripts/mkmpi.sh" \
  --compiler "$TASK_DIR/compiler dir/mpicc" --out "$output" 2>/dev/null; then
  echo 'mkmpi accepted a failed compiler' >&2; exit 1
fi
if [[ $(cat "$output") != 'previous provider' ]]; then
  echo 'mkmpi damaged the previous provider on failure' >&2; exit 1
fi
literal_flag='-DVALUE=$(never-run-this)'
"$ROOT/scripts/mkmpi.sh" --compiler "$TASK_DIR/compiler dir/mpicc" \
  --out "$output" --cflag "$literal_flag" --ldflag '-L/path with spaces'
if [[ $(cat "$output") != 'new provider' ]]; then
  echo 'mkmpi did not publish the provider' >&2; exit 1
fi
for expected in "$literal_flag" '-L/path with spaces' '-shared' \
  "-I$PLUMED_ROOT/include/plumed/wrapper"; do
  if ! grep -Fx -- "$expected" "$ARGUMENT_LOG" >/dev/null; then
    echo "mkmpi changed argument: $expected" >&2; exit 1
  fi
done
if [[ $(uname -s) == Darwin ]] && \
    ! grep -Fx -- '-Wl,-install_name,@rpath/provider.so' "$ARGUMENT_LOG" >/dev/null; then
  echo 'mkmpi used a temporary directory as the dylib install name' >&2; exit 1
fi
if "$ROOT/scripts/mkmpi.sh" --compiler relative/mpicc --out "$output" 2>/dev/null; then
  echo 'mkmpi accepted a relative compiler path' >&2; exit 1
fi
if "$ROOT/scripts/mkmpi.sh" --compiler 2>/dev/null; then
  echo 'mkmpi accepted a missing option argument' >&2; exit 1
fi
if grep -E -- '^-std=|^-x$' "$ARGUMENT_LOG" >/dev/null; then
  echo 'mkmpi added unconfigured language flags' >&2; exit 1
fi
shopt -s nullglob
temporary=("$TASK_DIR/output dir"/.plumed-mkmpi.*)
if ((${#temporary[@]})); then
  echo 'mkmpi left temporary build directories behind' >&2; exit 1
fi
echo 'mkmpi public-header arguments, cleanup, and publication passed'
