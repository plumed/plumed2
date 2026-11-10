#!/usr/bin/env python3
"""Benchmark IMAGE_DISTANCE against its explicitly selected serial Git baseline."""

import argparse
import csv
import math
import os
import shutil
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--reference-ref", required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--launcher", type=Path, required=True)
    parser.add_argument("--openmpi-local", action="store_true")
    args = parser.parse_args()
    checkout = Path(__file__).resolve().parents[3]
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    source = subprocess.check_output(
        ["git", "show", f"{args.reference_ref}:src/colvar/ImageDistance.cpp"],
        cwd=checkout,
        text=True,
    )
    source = source.replace('#include "Colvar.h"', '#include "plumed/core/Colvar.h"')
    source = source.replace('#include "core/', '#include "plumed/core/')
    source = source.replace('#include "tools/', '#include "plumed/tools/')
    source = source.replace("ImageDistance", "ImageDistanceReference")
    source = source.replace('"IMAGE_DISTANCE"', '"IMAGE_DISTANCE_REFERENCE"')
    build = work / "build"
    build.mkdir(exist_ok=True)
    (build / "reference.cpp").write_text(source)
    test = Path(__file__).resolve().parent
    shutil.copyfile(test / "main.cpp", build / "main.cpp")
    shutil.copyfile(
        test.parent / "rt-image-distance" / "test-tools.h", build / "test-tools.h"
    )
    shutil.copyfile(checkout / "regtest/scripts/mpi-test.h", build / "mpi-test.h")
    # Reuse exactly the compiler, flags and kernel objects of the built fixture.
    # Compile the reference into this executable too: a separately linked plugin
    # would register in a different kernel registry for static test executables.
    shutil.copyfile(args.exe.resolve().parent / "Makefile", build / "Makefile")
    includes = str(checkout / "src/include")
    if os.environ.get("CPLUS_INCLUDE_PATH"):
        includes += os.pathsep + os.environ["CPLUS_INCLUDE_PATH"]
    build_env = dict(os.environ, CPLUS_INCLUDE_PATH=includes)
    with (work / "reference-build.log").open("w") as log:
        subprocess.run(
            ["make", "exe", "CXX_SRC=main.cpp reference.cpp"],
            cwd=build,
            env=build_env,
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    executable = build / "exe"
    cases = [
        ("reference", 1, 1),
        ("serial", 1, 1),
        ("parallel", 1, 1),
        ("parallel", 1, 2),
        ("parallel", 1, 4),
        ("parallel", 1, 8),
        ("parallel", 2, 1),
        ("parallel", 4, 1),
        ("parallel", 2, 2),
        ("parallel", 4, 2),
    ]
    reference_values = {}
    with (work / "results.csv").open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(
            ["mode", "ranks", "threads", "atoms", "frames", "median_ms", "distance"]
        )
        for mode, ranks, threads in cases:
            run = work / f"{mode}-{ranks}x{threads}"
            run.mkdir(exist_ok=True)
            env = dict(
                os.environ,
                PLUMED_NUM_THREADS=str(threads),
                OMP_NUM_THREADS=str(threads),
                OMP_DYNAMIC="FALSE",
                PLUMED_MAXBACKUP="-1",
                IMAGE_DISTANCE_REFERENCE_BUILTIN="1",
            )
            command = [str(args.launcher.resolve())]
            if args.openmpi_local:
                command += ["--host", f"localhost:{ranks}", "--bind-to", "none"]
            command += [
                "-n",
                str(ranks),
                str(executable),
                "-n",
                str(ranks),
                "--benchmark",
                mode,
            ]
            completed = subprocess.run(
                command,
                cwd=run,
                env=env,
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            (run / "timing.log").write_text(completed.stdout)
            completed.check_returncode()
            rows = [
                line.split(",")
                for line in completed.stdout.splitlines()
                if line.startswith(f"{mode},")
            ]
            if len(rows) != 3:
                raise RuntimeError(
                    f"Expected three timing rows for {mode} {ranks}x{threads}"
                )
            for row in rows:
                atoms, distance = int(row[3]), float(row[6])
                if mode == "reference":
                    if not math.isfinite(distance):
                        raise RuntimeError("Nonfinite reference distance")
                    reference_values[atoms] = distance
                elif not math.isclose(
                    distance, reference_values[atoms], rel_tol=0, abs_tol=1e-12
                ):
                    raise RuntimeError(f"Distance changed in {mode} {ranks}x{threads}")
            writer.writerows(rows)
            output.flush()
            for row in rows:
                print(",".join(row), flush=True)


if __name__ == "__main__":
    main()
