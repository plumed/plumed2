#!/usr/bin/env bash
# Validate benchmark's actual MPI calculation, not its timing measurements.
set -euo pipefail
if (($# != 1)); then
  echo 'usage: check-thread-mpi-benchmark.sh ABSOLUTE_PLUMED_PROGRAM' >&2
  exit 2
fi
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-thread-benchmark.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
cp "$ROOT/regtest/basic/rt-benchmark/plumed.dat" "$TASK_DIR/plumed.dat"
python3 - "$1" "$TASK_DIR" <<'PYTHON'
import pathlib
import subprocess
import sys

program, directory = sys.argv[1:]
arguments = ["benchmark", "--plumed", "plumed.dat", "--natoms", "64",
             "--nsteps", "4", "--atom-distribution", "sc", "--domain-decomposition"]
expected = None
for ranks in (0, 1, 2, 4):
    launch = ["--no-mpi"] if ranks == 0 else ["--tmpi-threads", str(ranks)]
    result = subprocess.run([program] + launch + arguments, cwd=directory,
                            text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=30)
    if result.returncode:
        print(result.stdout, file=sys.stderr)
        raise SystemExit(f"benchmark failed with {ranks} ranks: {result.returncode}")
    path = pathlib.Path(directory, "colvar")
    actual = path.read_text()
    if expected is None:
        expected = actual
        rows = [line for line in actual.splitlines() if not line.startswith("#")]
        if len(rows) != 4:
            raise SystemExit(f"serial benchmark produced {len(rows)} rows instead of 4")
    elif actual != expected:
        raise SystemExit(f"benchmark collective results differ from serial with {ranks} ranks")
    path.unlink()
print("benchmark domain-decomposition results match serial at 1/2/4 thread-MPI ranks")
PYTHON
