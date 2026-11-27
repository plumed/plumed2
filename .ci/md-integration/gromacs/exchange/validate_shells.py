"""Check converged shell iterations freeze ABMD, including nested kernels."""

import argparse
from collections import defaultdict
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--plumed")
    parser.add_argument("--backend", choices=["mpi", "tmpi"], default="mpi")
    parser.add_argument("--launcher", default="")
    options = parser.parse_args()
    kernel = Path(options.kernel).resolve()
    gmx = str(Path(options.gmx).resolve())
    plumed = options.plumed or str(kernel.parent.parent / "bin/plumed")
    root = Path(options.output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    env = os.environ | {"PLUMED_KERNEL": str(kernel)}
    fixture = Path(__file__).resolve().parent / "shell"

    def run(command, cwd, name):
        with (cwd / (name + ".log")).open("w") as log:
            subprocess.run(
                command,
                cwd=cwd,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=120,
            )

    shutil.copy2(fixture / "ShellReviewTrace.cpp", root)
    run([plumed, "--no-mpi", "mklib", "ShellReviewTrace.cpp"], root, "mklib")
    library = root / (
        "ShellReviewTrace.dylib" if sys.platform == "darwin" else "ShellReviewTrace.so"
    )
    results = []
    for mode in ["flat", "local", "loaded"]:
        work = root / mode
        work.mkdir(exist_ok=True)
        for name in ["conf.gro", "md.mdp", "topol.top"]:
            shutil.copy2(fixture / name, work)
        guest = (
            f"LOAD FILE={library}\n"
            "d: DISTANCE ATOMS=1,2 NOPBC\n"
            "a: ABMD ARG=d TO=0 KAPPA=1000\n"
            "t: SHELLREVIEWTRACE ARG=d,a.d_min,a.bias\n"
            "b: BIASVALUE ARG=t\n"
        )
        (work / "guest.dat").write_text(guest)
        outer = guest if mode == "flat" else "p: PLUMED FILE=guest.dat"
        if mode == "loaded":
            outer += f" KERNEL={kernel}"
        (work / "plumed.dat").write_text(outer + "\n")
        run(
            [
                gmx,
                "grompp",
                "-f",
                "md.mdp",
                "-c",
                "conf.gro",
                "-p",
                "topol.top",
                "-o",
                "topol.tpr",
            ],
            work,
            "grompp",
        )
        launcher = shlex.split(options.launcher.format(n=1))
        command = launcher + [
            gmx,
            "mdrun",
            "-s",
            "topol.tpr",
            "-deffnm",
            "run",
            "-plumed",
            "plumed.dat",
            "-ntomp",
            "1",
            "-nb",
            "cpu",
            "-pme",
            "cpu",
            "-update",
            "cpu",
            "-pin",
            "off",
        ]
        if options.backend == "tmpi":
            command += ["-ntmpi", "1"]
        run(command, work, "mdrun")
        rows = defaultdict(list)
        pattern = r"SHELL_REVIEW step=(\d+) auxiliary=(\d) d=(\S+) min=(\S+) bias=(\S+)"
        for match in re.finditer(pattern, (work / "mdrun.log").read_text()):
            step, auxiliary, distance, minimum, bias = match.groups()
            values = [float(distance), float(minimum), float(bias)]
            if not all(math.isfinite(v) for v in values):
                raise AssertionError(f"{mode}: nonfinite shell evaluation")
            rows[int(step)].append((int(auxiliary), *values))
        if set(rows) != {0, 1, 2} or len(rows[0]) < 2:
            raise AssertionError(f"{mode}: missing converged shell evaluations")
        if abs(rows[0][0][1] - rows[0][-1][1]) < 0.05:
            raise AssertionError(f"{mode}: shell fixture did not change coordinates")
        for step, calls in rows.items():
            if calls[0][0] != 0 or any(c[0] != 1 for c in calls[1:]):
                raise AssertionError(f"{mode}: wrong auxiliary context at step {step}")
            if any(abs(c[2] - calls[0][2]) > 1e-12 for c in calls[1:]):
                raise AssertionError(
                    f"{mode}: ABMD minimum changed during shell iterations"
                )
        results.append({"case": mode, "status": "pass", "calls": dict(rows)})
        print(f"{mode}: converged shell iterations preserve ABMD state", flush=True)
    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
