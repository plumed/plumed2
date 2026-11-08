#!/usr/bin/env python3
"""Run normally, then exercise Abort dispatch in an isolated process.

Examples:
  python3 check-isolated.py ./exe
  python3 check-isolated.py mpiexec -n 4 ./exe
The same command is used for both runs; the executable receives -n 4 itself.
"""
import subprocess
import sys


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    arguments = sys.argv[1:]
    abort_only = arguments[0] == "--abort-only"
    if abort_only:
        arguments = arguments[1:]
    if not arguments:
        return 2
    command = arguments + ["-n", "4"]
    try:
        if not abort_only:
            normal = subprocess.run(command, timeout=120, check=False)
            if normal.returncode != 0:
                return 1
        abort = subprocess.run(command + ["--abort"], timeout=30, check=False,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except subprocess.TimeoutExpired:
        print("MPI API process did not finish within its timeout", file=sys.stderr)
        return 1
    # Runtime crash diagnostics need not be UTF-8. Match our ASCII markers
    # directly in the bytes; decode only for displaying failures.
    # A native Abort may return: MPI only promises a best-effort termination.
    # In that case the probe exits explicitly, rather than stranding its peers.
    # This checks dispatch and isolated completion, not guaranteed native abort.
    if (abort.returncode == 0 or b"MPI API intentional abort" not in abort.stdout or
            b"MPI_Abort unexpectedly threw" in abort.stdout):
        print(abort.stdout.decode("utf-8", errors="backslashreplace"), file=sys.stderr)
        print("MPI_Abort isolated dispatch probe failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
