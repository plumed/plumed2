#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiwrap-regtest.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/regtest/scripts" "$tmp/regtest/case" "$tmp/regtest/case/trajectories"
cp "$repo/regtest/scripts/run" "$tmp/regtest/scripts/run"
cp "$repo/regtest/scripts/fixzeroes.sh" "$tmp/regtest/scripts/fixzeroes.sh"
cp "$repo/regtest/scripts/exe.make" "$tmp/regtest/scripts/exe.make"
cp "$repo/regtest/scripts/mpi-cli.cpp" "$tmp/regtest/scripts/mpi-cli.cpp"
: > "$tmp/regtest/case/trajectories/testUtils.h"
mkdir -p "$tmp/bin"
cat > "$tmp/bin/plumed" <<'EOF'
#!/usr/bin/env bash
if [[ "$*" = "--no-mpi info --root" ]]; then
  echo "$TEST_REPO_ROOT"
  exit 0
fi
if [[ "$*" = "--no-mpi info --soext" ]]; then
  echo so
  exit 0
fi
if [[ "$*" = "--has-mpi" ]]; then
  printf 'QUERY %s\n' "$*" >> "$TEST_QUERY_LOG"
  if [[ ${TEST_NATIVE_MPI:-no} = yes || ${TEST_THREAD_MPI:-no} = yes ]]; then
    exit 0
  fi
  exit 1
fi
if [[ "$*" = "--no-mpi --has-real-mpi" ]]; then
  [[ ${TEST_INJECTED_REAL_MPI:-no} = yes ]]
  exit
fi
if [[ "$*" = "--no-mpi config -q has real_mpi" ]]; then
  if [[ ${TEST_OLD_CONFIG:-no} = yes ]]; then
    echo 'unknown feature real_mpi' >&2
    exit 2
  fi
  [[ ${TEST_NATIVE_MPI:-no} = yes ]]
  exit
fi
if [[ "$*" = "--no-mpi info --configuration" ]]; then
  if [[ ${TEST_NATIVE_MPI:-no} = yes ]]; then
    echo 'CPPFLAGS=-D__PLUMED_HAS_MPI'
  else
    echo 'CPPFLAGS=-D__PLUMED_MPIWRAP'
  fi
  exit 0
fi
if [[ "$*" = "--no-mpi config -q mpiexec" || "$*" = "--no-mpi config mpiexec" ]]; then
  [[ ${TEST_NATIVE_MPI:-no} = yes ]] || exit 1
  echo "$TEST_CONFIG_MPIRUN"
  exit 0
fi
if [[ "$*" = "--has-thread-mpi" || "$*" = "--no-mpi --has-thread-mpi" ]]; then
  printf 'QUERY %s\n' "$*" >> "$TEST_QUERY_LOG"
  if [[ ${TEST_OLD_CONFIG:-no} = yes ]]; then
    echo 'unknown option --has-thread-mpi' >&2
    exit 2
  fi
  [[ ${TEST_THREAD_MPI:-no} = yes ]]
  exit
fi
if [[ "$*" = "--is-installed" || "$*" = "--no-mpi --is-installed" ]]; then
  [[ ${TEST_INSTALLED:-no} = yes ]]
  exit
fi
if [[ " $* " = *" driver "* || " $* " = *" sum_hills "* || " $* " = *" benchmark "* || " $* " = *" simplemd "* ]]; then
  printf 'PLUMED %s\n' "$*" >> "$TEST_LAUNCH_LOG"
  exit 0
fi
echo "unexpected fake plumed arguments: $*" >&2
exit 2
EOF
cat > "$tmp/bin/fake-mpiexec" <<'EOF'
#!/usr/bin/env bash
printf 'MPI %s\n' "$*" >> "$TEST_LAUNCH_LOG"
[[ ${1:-} = -np && $# -ge 4 ]] || exit 2
shift 2
exec "$@"
EOF
cat > "$tmp/bin/config-mpiexec" <<'EOF'
#!/usr/bin/env bash
printf 'CONFIG %s\n' "$*" >> "$TEST_LAUNCH_LOG"
[[ ${1:-} = -np && $# -ge 4 ]] || exit 2
shift 2
exec "$@"
EOF
cat > "$tmp/bin/fake-make" <<'EOF'
#!/usr/bin/env bash
printf 'MAKE %s\n' "$*" >> "$TEST_LAUNCH_LOG"
if [[ " $* " = *" exe "* ]]; then
  cat > exe <<'EXE'
#!/usr/bin/env bash
printf 'MPI_EXE %s\n' "$*" >> "$TEST_LAUNCH_LOG"
EXE
  chmod +x exe
fi
EOF
chmod +x "$tmp/bin/plumed" "$tmp/bin/fake-mpiexec" "$tmp/bin/config-mpiexec" "$tmp/bin/fake-make"

run_case() {
  local name=$1 mpi=$2 threads=$3 per_test=$4 chdir=$5 override=$6 expect=$7 old=$8
  local test_type=${9:-driver}
  local program_name=${10:-plumed}
  local cli_program_name=${11:-}
  local test_arguments='--test-ranks 3'
  if [[ ${12:-no} = yes ]]; then test_arguments=''; fi
  if [[ "$test_type" = plumed ]]; then test_arguments="benchmark $test_arguments"; fi
  local dir="$tmp/regtest/case/$name"
  mkdir -p "$dir"
  cat > "$dir/config" <<EOF
type=$test_type
mpiprocs=3
arg="$test_arguments"
plumed_thread_mpi=$per_test
plumed_chdir=$chdir
plumed_mpi_harness=${12:-no}
EOF
  if [[ "$test_type" = simplemd ]]; then echo 'nstep 1' > "$dir/in"; fi
  if [[ "$test_type" = make ]]; then echo 'int main() { return 0; }' > "$dir/test.cpp"; fi
  : > "$TEST_LAUNCH_LOG"
  : > "$TEST_QUERY_LOG"
  if ! (cd "$dir" && TEST_NATIVE_MPI="$mpi" TEST_THREAD_MPI="$threads" TEST_OLD_CONFIG="$old" \
    TEST_REPO_ROOT="$repo" TEST_CONFIG_MPIRUN="$tmp/bin/config-mpiexec" \
    TEST_LAUNCH_LOG="$TEST_LAUNCH_LOG" PATH="$tmp/bin:$PATH" \
    PLUMED_PREPEND_PATH="$tmp/bin" PLUMED_MPIRUN="$override" PLUMED_MAKE="$tmp/bin/fake-make" \
    PLUMED_PROGRAM_NAME="$program_name" PLUMED_CLI_PROGRAM_NAME="$cli_program_name" \
    "$tmp/regtest/scripts/run" >run.out 2>run.err); then
    cat "$dir/run.out" "$dir/run.err" >&2
    return 1
  fi
  if [[ -n $expect ]]; then
    if ! grep -q "$expect" "$dir/report.txt"; then
      cat "$dir/report.txt" >&2
      return 1
    fi
  fi
}

check_log() {
  local pattern=$1
  if ! grep -q "$pattern" "$TEST_LAUNCH_LOG"; then
    cat "$TEST_LAUNCH_LOG" >&2
    return 1
  fi
}

# A CLI-only injection host advertises real MPI, unlike the serial executable.
cat > "$tmp/bin/mpi-host" <<'EOF'
#!/usr/bin/env bash
printf 'HOST %s\n' "$*" >> "$TEST_LAUNCH_LOG"
TEST_INJECTED_REAL_MPI=yes TEST_NATIVE_MPI=no TEST_THREAD_MPI=yes exec "$(dirname "$0")/plumed" "$@"
EOF
chmod +x "$tmp/bin/mpi-host"

export TEST_LAUNCH_LOG="$tmp/launch.log"
export TEST_QUERY_LOG="$tmp/query.log"
run_case internal no yes auto no "" "Using PLUMED internal thread-MPI with 3 ranks" no
check_log 'PLUMED --tmpi-threads 3 driver --test-ranks 3'
grep -q '^QUERY --has-mpi$' "$TEST_QUERY_LOG"
grep -q '^QUERY --no-mpi --has-thread-mpi$' "$TEST_QUERY_LOG"

run_case benchmark no yes auto no "" "Using PLUMED internal thread-MPI with 3 ranks" no plumed
check_log 'PLUMED --tmpi-threads 3 benchmark --test-ranks 3'

run_case sum_hills no yes auto no "" "Using PLUMED internal thread-MPI with 3 ranks" no sum_hills
check_log 'PLUMED --tmpi-threads 3 sum_hills --test-ranks 3'

run_case simplemd no yes auto no "" "Using PLUMED internal thread-MPI with 3 ranks" no simplemd
check_log 'PLUMED --tmpi-threads 3 simplemd in'

run_case external no yes auto no "" "MPI NOT INSTALLED" no make
test ! -s "$TEST_LAUNCH_LOG"

run_case make-internal no yes auto no "" "Using PLUMED internal thread-MPI with 3 ranks" no make plumed "" yes
check_log '^MPI_EXE -n 3$'

run_case make-native yes no auto no "$tmp/bin/fake-mpiexec" "" no make plumed "" yes
check_log '^MPI -np 3 env ./exe -n 3$'

run_case native yes no auto no "" "" no
check_log '^CONFIG -np 3 '
check_log '^PLUMED driver --test-ranks 3$'
grep -q '^QUERY --has-mpi$' "$TEST_QUERY_LOG"

run_case old_native yes no auto no "" "" yes
check_log '^CONFIG -np 3 '
check_log '^PLUMED driver --test-ranks 3$'
grep -q '^QUERY --has-mpi$' "$TEST_QUERY_LOG"
grep -q '^QUERY --no-mpi --has-thread-mpi$' "$TEST_QUERY_LOG"

run_case native_override yes no auto no "$tmp/bin/fake-mpiexec" "" no
check_log '^MPI -np 3 '
if grep -q '^CONFIG ' "$TEST_LAUNCH_LOG"; then
  echo "explicit PLUMED_MPIRUN was ignored" >&2
  exit 1
fi

run_case native_absolute_program yes no auto no "" "" no driver "$tmp/bin/plumed"
check_log '^CONFIG -np 3 '
check_log '^PLUMED driver --test-ranks 3$'

run_case old_serial no no auto no "" "MPI NOT INSTALLED" yes
test ! -s "$TEST_LAUNCH_LOG"
grep -q '^QUERY --no-mpi --has-thread-mpi$' "$TEST_QUERY_LOG"

run_case thread_disabled no yes no no "" "thread-MPI disabled by test config" no
test ! -s "$TEST_LAUNCH_LOG"

run_case thread_override no yes auto no "$tmp/bin/fake-mpiexec" "MPI NOT INSTALLED" no
test ! -s "$TEST_LAUNCH_LOG"

run_case chdir no yes no yes "" "CHDIR is unsupported with thread-MPI" no
test ! -s "$TEST_LAUNCH_LOG"

echo "regtest runner selection passed (CLI types, external executable, native, legacy, absolute program, thread fallback, override, CHDIR skip)"

# CLI commands use the injected host and native launcher; external-engine tests
# retain the serial executable and cannot accidentally borrow its injection.
run_case injected_cli no yes auto no "$tmp/bin/fake-mpiexec" "" no driver plumed "$tmp/bin/mpi-host"
check_log '^HOST driver --test-ranks 3$'
check_log '^MPI -np 3 '
grep -q '^HOST --no-mpi --has-real-mpi$' "$TEST_LAUNCH_LOG"
run_case injected_external no yes auto no "$tmp/bin/fake-mpiexec" "MPI NOT INSTALLED" no make plumed "$tmp/bin/mpi-host"
test ! -s "$TEST_LAUNCH_LOG"
echo 'CLI-only injection host selection checks passed'
