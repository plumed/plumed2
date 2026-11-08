#!/usr/bin/env bash
# Exercise argument boundaries, installation independence, and failed publication.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-mkmpi-test.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
TASK_DIR=$(cd "$TASK_DIR" && pwd)
export PLUMED_ROOT="$TASK_DIR/installed tree"
mkdir -p "$PLUMED_ROOT/src/config" "$PLUMED_ROOT/src/main" "$PLUMED_ROOT/include/plumed/wrapper" "$TASK_DIR/compiler dir" "$TASK_DIR/output dir"
printf '/* public header */\n' > "$PLUMED_ROOT/include/plumed/wrapper/Plumed.h"
export PLUMED_INCLUDEDIR="$PLUMED_ROOT/include"
export PLUMED_PROGRAM_NAME=plumed
printf "mpi_provider_compile_flags='-fPIC'\nmpi_provider_link_flags='-shared'\nsoext=so\n" > "$PLUMED_ROOT/src/config/compile_options.sh"
cp "$ROOT/src/main/main.cpp" "$PLUMED_ROOT/src/main/main.cpp"
printf 'kernel\n' > "$TASK_DIR/libplumedKernel.so"
printf 'wrapper\n' > "$TASK_DIR/libplumed.so"
export ARGUMENT_LOG="$TASK_DIR/arguments" COMPILER_LOG="$TASK_DIR/compiler-used"
cat > "$TASK_DIR/compiler dir/mpicc" <<'COMPILER'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$@" > "$ARGUMENT_LOG"
printf '%s\n' "$0" > "$COMPILER_LOG"
output=
while (($#)); do
  if [[ "${1##*/}" = main.c ]]; then
    cmp "$1" "$PLUMED_ROOT/src/main/main.cpp"
  fi
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
for mode in shared runtime; do
  "$ROOT/scripts/mkmpi.sh" --executable --link "$mode" \
    --compiler "$TASK_DIR/compiler dir/mpicc" --out "$output"
  grep -E -- '/main\.c$' "$ARGUMENT_LOG" >/dev/null
  grep -Fx -- '-D__PLUMED_HAS_MPI=1' "$ARGUMENT_LOG" >/dev/null
  if grep -E -- '^-std=|^-x$' "$ARGUMENT_LOG" >/dev/null; then
    echo 'mkmpi forced an executable language mode' >&2; exit 1
  fi
  if grep -Fx -- '-shared' "$ARGUMENT_LOG" >/dev/null; then
    echo 'mkmpi linked its executable as a provider library' >&2; exit 1
  fi
  if [[ "$mode" = shared ]]; then
    grep -Fx -- '-lplumed' "$ARGUMENT_LOG" >/dev/null
  else
    grep -Fx -- '-D__PLUMED_WRAPPER_LINK_RUNTIME=1' "$ARGUMENT_LOG" >/dev/null
    if grep -Fx -- '-lplumed' "$ARGUMENT_LOG" >/dev/null; then
      echo 'runtime executable unexpectedly links libplumed' >&2; exit 1
    fi
  fi
done
for invalid in runtime other; do
  if "$ROOT/scripts/mkmpi.sh" --link "$invalid" --out "$output" 2>/dev/null; then
    echo 'mkmpi accepted a link mode without executable mode' >&2; exit 1
  fi
done
# Isolate discovery from MPI installations on the developer's PATH.
mkdir "$TASK_DIR/bin"
for program in bash basename cat cmp cp dirname mktemp mv rm uname; do
  ln -s "$(command -v "$program")" "$TASK_DIR/bin/$program"
done
for name in mpicc mpicxx env-cxx; do
  cp "$TASK_DIR/compiler dir/mpicc" "$TASK_DIR/bin/$name"
done
check_default() {
  local expected=$1
  shift
  env PATH="$TASK_DIR/bin" MPICC= MPICXX= "$@" \
    "$ROOT/scripts/mkmpi.sh" --out "$output"
  if [[ $(cat "$COMPILER_LOG") != "$expected" ]]; then
    echo "mkmpi selected the wrong default compiler (expected $expected)" >&2
    exit 1
  fi
}
check_default "$TASK_DIR/compiler dir/mpicc" \
  MPICC="$TASK_DIR/compiler dir/mpicc" MPICXX=env-cxx
(cd "$TASK_DIR" && check_default "$TASK_DIR/compiler dir/mpicc" MPICC="compiler dir/mpicc")
check_default "$TASK_DIR/bin/env-cxx" MPICC=/missing/mpicc MPICXX=env-cxx
check_default "$TASK_DIR/bin/mpicc"
# Executable mode uses the same C-first compiler preference as provider mode.
env PATH="$TASK_DIR/bin" MPICC=mpicc MPICXX=env-cxx \
  "$ROOT/scripts/mkmpi.sh" --executable --out "$output"
test "$(cat "$COMPILER_LOG")" = "$TASK_DIR/bin/mpicc"
rm "$TASK_DIR/bin/mpicc"
check_default "$TASK_DIR/bin/mpicxx"
# Explicit selection must take precedence over all defaults.
MPICC=/missing/mpicc MPICXX=env-cxx PATH="$TASK_DIR/bin" \
  "$ROOT/scripts/mkmpi.sh" --compiler "$TASK_DIR/compiler dir/mpicc" --out "$output"
if [[ $(cat "$COMPILER_LOG") != "$TASK_DIR/compiler dir/mpicc" ]]; then
  echo 'mkmpi ignored the explicit compiler' >&2; exit 1
fi
rm "$TASK_DIR/bin/mpicxx"
if env PATH="$TASK_DIR/bin" MPICC= MPICXX= \
    "$ROOT/scripts/mkmpi.sh" --out "$output" >"$TASK_DIR/no-compiler.log" 2>&1; then
  echo 'mkmpi accepted missing compiler defaults' >&2; exit 1
fi
grep -q 'no MPI compiler found' "$TASK_DIR/no-compiler.log"
shopt -s nullglob
temporary=("$TASK_DIR/output dir"/.plumed-mkmpi.*)
if ((${#temporary[@]})); then
  echo 'mkmpi left temporary build directories behind' >&2; exit 1
fi
echo 'mkmpi compiler discovery, public-header arguments, cleanup, and publication passed'
