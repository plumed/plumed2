#!/usr/bin/env python3
"""Exercise AmberTools sander against a serial PLUMED kernel."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

KAPPA = 1000.0
AT = 0.2


def read_colvar(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        row = [float(value) for value in line.split()]
        if len(row) != 4 or not all(math.isfinite(value) for value in row):
            raise RuntimeError(f"Malformed or non-finite COLVAR row in {path}: {line}")
        rows.append(row)
    if not rows:
        raise RuntimeError(f"No data rows in {path}")
    return rows


def read_mdcrd(path: Path, natoms: int = 12) -> list[list[float]]:
    """Read the 8-column ASCII Amber trajectory written by ioutfm=0."""
    lines = path.read_text().splitlines()[1:]
    values: list[float] = []
    for line in lines:
        values.extend(
            float(line[i : i + 8])
            for i in range(0, len(line), 8)
            if line[i : i + 8].strip()
        )
    frame_size = 3 * natoms
    if (
        not values
        or len(values) % frame_size
        or not all(math.isfinite(value) for value in values)
    ):
        raise RuntimeError(f"Unexpected coordinate count in {path}: {len(values)}")
    return [values[i : i + frame_size] for i in range(0, len(values), frame_size)]


def read_restart(path: Path, natoms: int = 12) -> list[float]:
    """Read the coordinate block from an ASCII Amber restart file."""
    lines = path.read_text().splitlines()[2:]
    values: list[float] = []
    for line in lines:
        values.extend(
            float(line[i : i + 12])
            for i in range(0, len(line), 12)
            if line[i : i + 12].strip()
        )
        if len(values) >= 3 * natoms:
            return values[: 3 * natoms]
    raise RuntimeError(f"Could not read {3 * natoms} coordinates from {path}")


def run_case(
    *, root: Path, label: str, nprocs: int, enabled: bool, args: argparse.Namespace
) -> dict:
    case = root / label
    case.mkdir()
    for name in ("NMA.prmtop", "NMA.inpcrd"):
        shutil.copy2(args.topology_dir / name, case / name)
    shutil.copy2(
        args.fixture_dir / ("mdin" if enabled else "mdin-no-plumed"), case / "mdin"
    )
    if enabled:
        shutil.copy2(args.fixture_dir / "plumed.dat", case / "plumed.dat")

    if args.mode == "serial":
        command = []
    else:
        command = [args.launcher]
        if args.launcher_kind == "openmpi":
            command += ["--oversubscribe", "--bind-to", "none", "-n", str(nprocs)]
        else:
            command += ["-n", str(nprocs)]
    command += [
        str(args.sander),
        "-O",
        "-i",
        "mdin",
        "-p",
        "NMA.prmtop",
        "-c",
        "NMA.inpcrd",
        "-o",
        "mdout",
        "-r",
        "restart.rst",
        "-x",
        "mdcrd",
        "-frc",
        "mdfrc",
    ]
    env = os.environ.copy()
    env.update(
        OMP_NUM_THREADS="1",
        PLUMED_NUM_THREADS="1",
        PLUMED_KERNEL=str(args.kernel)
        if enabled
        else "/definitely/not/a/plumed-kernel",
        MPICC="/usr/bin/false",
        MPICXX="/usr/bin/false",
    )
    if args.mode == "mpi":
        env.update(
            PLUMED_MPI_WRAPPER="/definitely/not/a/plumed-mpi-provider",
            PLUMED_MPI_WRAPPER_FORCE="0",
        )
    else:
        env.pop("PLUMED_MPI_WRAPPER", None)
        env.pop("PLUMED_MPI_WRAPPER_FORCE", None)
    log = subprocess.run(
        command,
        cwd=case,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (case / "command.txt").write_text(
        shlex.join([os.fspath(value) for value in command]) + "\n"
    )
    (case / "run.log").write_text(log.stdout)
    if log.returncode:
        raise RuntimeError(
            f"{label} exited with {log.returncode}; see {case / 'run.log'}\n"
            + "\n".join(log.stdout.splitlines()[-60:])
        )

    mdout = (case / "mdout").read_text(errors="replace")
    if args.mode == "mpi":
        amber_ranks = re.findall(
            r"Running AMBER/MPI version on\s+(\d+)\s+nodes?\b", mdout
        )
        if amber_ranks != [str(nprocs)]:
            raise RuntimeError(
                f"{label}: Amber reported native rank counts {amber_ranks}, expected exactly {nprocs}"
            )
    elif "Running AMBER/MPI version" in mdout:
        raise RuntimeError(f"{label}: serial sander printed the MPI engine banner")
    coords = read_mdcrd(case / "mdcrd")
    forces = read_mdcrd(case / "mdfrc")
    final = read_restart(case / "restart.rst")
    if len(coords) != 5 or len(forces) != 5 or not all(math.isfinite(x) for x in final):
        raise RuntimeError(f"{label}: malformed trajectory, forces, or restart")
    if enabled:
        if log.stdout.count("PLUMED: Finished setup") != 1:
            raise RuntimeError(
                f"{label} did not initialize PLUMED; see {case / 'run.log'}"
            )
        plumed_ranks = re.findall(
            r"PLUMED: Running over\s+(\d+)\s+nodes?\b", log.stdout
        )
        if plumed_ranks != [str(nprocs)]:
            raise RuntimeError(
                f"{label} did not report {nprocs} PLUMED ranks; see {case / 'run.log'}"
            )
        if not (case / "COLVAR").exists():
            raise RuntimeError(f"{label} did not produce COLVAR")
        rows = read_colvar(case / "COLVAR")
        if len(rows) != 6:
            raise RuntimeError(
                f"{label}: expected initial frame plus five COLVAR steps"
            )
        for step, row in enumerate(rows):
            time, distance, bias, force2 = row
            # Sander evaluates the initial forces once during setup, then
            # again at the first MD step before advancing its step counter.
            if abs(time - max(0, step - 1) * 0.001) > 1.0e-8:
                raise RuntimeError(f"{label}: missing or unordered COLVAR time step")
            expected_bias = 0.5 * KAPPA * (distance - AT) ** 2
            if abs(bias - expected_bias) > 1.0e-4:
                raise RuntimeError(
                    f"{label}: bias {bias} disagrees with analytic {expected_bias}"
                )
            if abs(force2 - 2.0 * KAPPA * bias) > 0.002:
                raise RuntimeError(
                    f"{label}: restraint force magnitude is inconsistent with its bias"
                )
        initial = read_restart(args.topology_dir / "NMA.inpcrd")
        for row, coordinates in zip(rows, [initial, initial] + coords[:-1]):
            native_distance = (
                math.sqrt(
                    sum((coordinates[i] - coordinates[i + 3]) ** 2 for i in range(3))
                )
                / 10.0
            )
            # Native ASCII coordinates have three decimals in angstrom.
            if abs(row[1] - native_distance) > 0.0002:
                raise RuntimeError(
                    f"{label}: COLVAR distance disagrees with Amber coordinates"
                )
        if len(coords) != 5 or len(final) != 36:
            raise RuntimeError(
                f"{label}: expected 5 trajectory frames and 36 final coordinates"
            )
        return {
            "label": label,
            "nprocs": nprocs,
            "colvar": rows,
            "trajectory": coords,
            "forces": forces,
            "final": final,
        }

    if (
        "PLUMED is on" in mdout
        or "PLUMED: " in log.stdout
        or (case / "COLVAR").exists()
    ):
        raise RuntimeError(f"{label} activated PLUMED while plumed=0")
    return {
        "label": label,
        "nprocs": nprocs,
        "trajectory": coords,
        "forces": forces,
        "final": final,
    }


def max_delta(a: list[float], b: list[float]) -> float:
    if len(a) != len(b):
        raise RuntimeError(f"Cannot compare vectors with lengths {len(a)} and {len(b)}")
    return max(abs(x - y) for x, y in zip(a, b))


def check_force_increment(active: dict, baseline: dict, initial: list[float]) -> float:
    """Check native Amber forces at identical initial coordinates, in kcal/mol/A."""
    separation = [initial[i] - initial[i + 3] for i in range(3)]
    length = math.sqrt(sum(x * x for x in separation))
    expected = [0.0] * 36
    # PLUMED defaults: kJ/mol and nm; Amber: kcal/mol and angstrom.
    for axis in range(3):
        expected[axis] = (
            -KAPPA * (length / 10.0 - AT) * separation[axis] / length / 41.84
        )
        expected[axis + 3] = -expected[axis]
    increment = [x - y for x, y in zip(active["forces"][0], baseline["forces"][0])]
    error = max_delta(increment, expected)
    # Two native force dumps rounded to three decimals each.
    if error > 0.0011:
        raise RuntimeError(f"Native initial force increment error: {error} kcal/mol/A")
    return error


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("mpi", "serial"), default="mpi")
    parser.add_argument("--launcher-kind", choices=("openmpi", "mpich"))
    parser.add_argument("--launcher", type=Path)
    parser.add_argument(
        "--sander", type=Path, required=True, help="Built sander executable"
    )
    parser.add_argument(
        "--kernel", type=Path, required=True, help="Serial PLUMED libplumedKernel"
    )
    parser.add_argument(
        "--results", type=Path, required=True, help="New output directory"
    )
    parser.add_argument(
        "--amber-source",
        type=Path,
        required=True,
        help="Extracted AmberTools source tree",
    )
    parser.add_argument(
        "--fixture-dir",
        type=Path,
        default=Path(__file__).resolve().with_name("input"),
        help="Directory containing mdin, mdin-no-plumed and plumed.dat",
    )
    parser.add_argument(
        "--ranks",
        type=int,
        nargs="+",
        choices=(1, 2, 4),
        default=[1, 2, 4],
        help="Rank counts to test (default: 1 2 4; use 1 for singleton smoke testing)",
    )
    args = parser.parse_args()
    if args.mode == "mpi" and (args.launcher is None or args.launcher_kind is None):
        parser.error("MPI mode requires --launcher and --launcher-kind")
    args.topology_dir = (
        args.amber_source
        / "AmberTools/test/paramfit/simplex_perfect_fit_NMA/prmtop_creation"
    )
    args.ranks = [1] if args.mode == "serial" else sorted(set(args.ranks))
    paths = [args.sander, args.kernel]
    if args.launcher is not None:
        paths.append(args.launcher)
    for path in paths:
        if not path.exists():
            parser.error(f"Path does not exist: {path}")
    if args.results.exists():
        parser.error(f"Results directory already exists: {args.results}")
    args.results.mkdir(parents=True)

    if 1 not in args.ranks:
        parser.error("--ranks must include 1 for the reference run")
    active = [
        run_case(
            root=args.results, label=f"active-{n}", nprocs=n, enabled=True, args=args
        )
        for n in args.ranks
    ]
    baseline = run_case(
        root=args.results, label="baseline-1", nprocs=1, enabled=False, args=args
    )
    if max(args.ranks) > 1:
        baseline_multi = run_case(
            root=args.results,
            label=f"baseline-{max(args.ranks)}",
            nprocs=max(args.ranks),
            enabled=False,
            args=args,
        )
        if len(baseline_multi["trajectory"]) != len(baseline["trajectory"]):
            raise RuntimeError(
                "Baseline trajectory frame counts differ between 1 and multi-rank runs"
            )
        for frame1, frameN in zip(baseline["trajectory"], baseline_multi["trajectory"]):
            if max_delta(frame1, frameN) > 0.0011:
                raise RuntimeError(
                    "Baseline trajectories differ between 1 rank and multi-rank runs"
                )

    reference = active[0]
    initial = read_restart(args.topology_dir / "NMA.inpcrd")
    force_error = check_force_increment(reference, baseline, initial)
    for case in active[1:]:
        if len(case["trajectory"]) != len(reference["trajectory"]):
            raise RuntimeError(
                f"{case['label']}: trajectory frame count differs from 1 rank"
            )
        for ref_frame, frame in zip(reference["trajectory"], case["trajectory"]):
            if max_delta(ref_frame, frame) > 0.005:
                raise RuntimeError(
                    f"{case['label']}: trajectory coordinates differ from 1-rank run"
                )
        check_force_increment(case, baseline, initial)
        for ref_force, force in zip(reference["forces"], case["forces"]):
            if max_delta(ref_force, force) > 0.002:
                raise RuntimeError(f"{case['label']}: forces differ from 1-rank run")
        if len(case["colvar"]) != len(reference["colvar"]):
            raise RuntimeError(
                f"{case['label']}: COLVAR frame count differs from 1 rank"
            )
        for ref_row, row in zip(reference["colvar"], case["colvar"]):
            if max_delta(ref_row, row) > 2.0e-5:
                raise RuntimeError(f"{case['label']}: COLVAR differs from 1-rank run")
        if max_delta(reference["final"], case["final"]) > 0.005:
            raise RuntimeError(
                f"{case['label']}: final coordinates differ from 1-rank run"
            )

    displacement = max_delta(reference["final"], baseline["final"])
    if displacement < 0.001:
        raise RuntimeError("PLUMED restraint did not change the final coordinates")
    summary = {
        "mode": args.mode,
        "frames_per_run": len(reference["colvar"]),
        "max_biased_vs_baseline_displacement_angstrom": displacement,
        "initial_distance_nm": reference["colvar"][0][1],
        "initial_bias_kj_per_mol": reference["colvar"][0][2],
        "initial_force_squared": reference["colvar"][0][3],
        "max_initial_atomic_force_error_kcal_per_mol_angstrom": force_error,
    }
    if args.mode == "mpi":
        summary["active_ranks"] = args.ranks
    else:
        summary["execution"] = "serial"
    (args.results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
