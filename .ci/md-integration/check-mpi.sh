#!/usr/bin/env bash
# Check the host MPI without loading either PLUMED or an MD engine.
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
"${MD_MPICXX:?Set MD_MPICXX}" "$root/regtest/basic/rt-mpiwrap-api/mpi-native-probe.cpp" \
  -o "$RUNNER_TEMP/md-native-probe"
python3 - <<'PY'
import os
import subprocess

command = [os.environ["MD_MPIEXEC"]]
if os.environ["MD_LAUNCHER_KIND"] == "openmpi":
    command += ["--oversubscribe", "--bind-to", "none"]
command += ["-n", "4", os.environ["RUNNER_TEMP"] + "/md-native-probe", "4"]
subprocess.run(command, check=True, timeout=60)
PY
