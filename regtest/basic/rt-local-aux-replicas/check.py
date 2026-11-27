"""Verify coupled auxiliary calls terminate with diagnostics instead of hanging."""

import os
from pathlib import Path
import shlex
import subprocess

launcher = (
    []
    if os.environ.get("PLUMED_TEST_INTERNAL") == "yes"
    else shlex.split(os.environ["PLUMED_TEST_MPI"])
)
for mode in (1, 2):
    command = launcher + ["./exe", "-n", "4", str(mode)]
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60
    )
    # Abort handlers and MPI runtimes may emit bytes that are not UTF-8.
    # Preserve the complete output, and check the required ASCII diagnostic
    # without decoding unrelated runtime/backtrace output.
    Path(f"negative-mode{mode}.log").write_bytes(result.stdout)
    if (
        result.returncode == 0
        or b"action ensemble (ENSEMBLE) communicates between replicas"
        not in result.stdout
    ):
        raise RuntimeError(f"mode {mode}: missing coupled-action diagnostic")
    with Path("output").open("a") as stream:
        stream.write(
            f"Coupled local auxiliary mode {mode} rejected before unmatched collectives: PASS\n"
        )
