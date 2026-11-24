#!/usr/bin/env bash
# Exercise unchanged communicator forwarding into a separately loaded kernel.
set -euo pipefail
if (($# != 2)); then
  echo 'usage: check-mpi-nested-kernel.sh ABSOLUTE_RUNTIME_PROGRAM ABSOLUTE_KERNEL' >&2
  exit 2
fi
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-nested-kernel.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
ulimit -c 0
cp "$2" "$TASK_DIR/guest-kernel.so"
cp "$ROOT/regtest/basic/rt-plumed-mpi/plumed.dat" "$TASK_DIR/plumed.dat"
cp "$ROOT/regtest/trajectories/trajectory.xyz" "$TASK_DIR/trajectory.xyz"
printf 'PLUMED FILE=plumed.dat KERNEL=%s\n' "$TASK_DIR/guest-kernel.so" > "$TASK_DIR/plumedx.dat"
python3 - "$1" "$2" "$TASK_DIR" "$ROOT" <<'PYTHON'
import os
import pathlib
import shutil
import subprocess
import sys

program, kernel, directory, root = sys.argv[1:]
environment = dict(os.environ, PLUMED_KERNEL=kernel)
for ranks in (1, 2, 4):
    result = subprocess.run(
        [program, "--tmpi-threads", str(ranks), "driver", "--plumed=plumedx.dat",
         "--timestep=0.05", "--ixyz", "trajectory.xyz", "--dump-forces", "ff",
         "--debug-pd", "yes"],
        cwd=directory, env=environment, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=30,
    )
    if result.returncode:
        print(result.stdout, file=sys.stderr)
        raise SystemExit(f"nested kernel failed with {ranks} ranks: {result.returncode}")
    actual = pathlib.Path(directory, "COLVAR").read_text()
    expected = pathlib.Path(root, "regtest/basic/rt-plumed-mpi/COLVAR.reference").read_text()
    if actual != expected:
        raise SystemExit(f"nested kernel COLVAR differs with {ranks} ranks")
    pathlib.Path(directory, "COLVAR").unlink()
# A single-replica PLUMED action runs its guest locally and does not forward an
# MPI communicator. Use the existing multi-replica fixture to exercise the
# actual GREX communicator-pointer boundary as well.
fixture = pathlib.Path(root, "regtest/basic/rt-plumed-mpi2")
for name in ("plumed.dat", "plumed.0.dat", "test.pdb"):
    shutil.copyfile(fixture / name, pathlib.Path(directory, name))
for replica in range(3):
    name = f"trajectory.{replica}.xyz"
    shutil.copyfile(pathlib.Path(root, "regtest/basic/rt-mpi3", name), pathlib.Path(directory, name))
result = subprocess.run(
    [program, "--tmpi-threads", "6", "driver", "--plumed=plumedx.dat",
     "--timestep=0.05", "--ixyz", "trajectory.xyz", "--dump-forces", "ff",
     "--pdb", "test.pdb", "--multi", "3", "--debug-dd", "yes",
     "--debug-grex", "2", "--debug-grex-log", "log.grex"],
    cwd=directory, env=environment, text=True, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT, timeout=30,
)
if result.returncode:
    print(result.stdout, file=sys.stderr)
    raise SystemExit(f"nested kernel replica communicator forwarding failed: {result.returncode}")
for replica in range(3):
    name = f"COLVAR.{replica}"
    actual = pathlib.Path(directory, name).read_text()
    expected = (fixture / (name + ".reference")).read_text()
    if actual != expected:
        raise SystemExit(f"nested kernel replica {replica} COLVAR differs")
print("nested guest kernel and replica communicator forwarding passed")
PYTHON
