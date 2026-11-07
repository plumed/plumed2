#!/usr/bin/env python3
"""Validate LAMMPS injection into serial PLUMED with real MPI and OpenMP."""

import argparse
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
MPI_CASES = [(1, 1), (2, 1), (4, 1), (4, 2)]
SERIAL_CASES = [(1, 1), (1, 2)]
NATOMS = 256
SPACING = 5.26
NSTEPS = 20
DT = 0.0005
TARGET = 3.5
KAPPA = 0.1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def write_data(path):
    """Assign IDs before decomposition; create_atoms IDs depend on rank count."""
    basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    with path.open("w") as out:
        out.write(f"Argon injection test\n\n{NATOMS} atoms\n1 atom types\n\n")
        for axis in "xyz":
            out.write(f"0 {4 * SPACING} {axis}lo {axis}hi\n")
        out.write("\nMasses\n\n1 39.948\n\nAtoms # atomic\n\n")
        atom = 0
        for cell in itertools.product(range(4), repeat=3):
            for offset in basis:
                atom += 1
                xyz = " ".join(
                    f"{SPACING * (i + j):.12f}" for i, j in zip(cell, offset)
                )
                out.write(f"{atom} 1 {xyz}\n")
        assert atom == NATOMS


def read_table(path):
    rows = [
        [float(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    require(all(math.isfinite(x) for row in rows for x in row), f"nonfinite {path}")
    return rows


def read_dump(path):
    frames = []
    lines = iter(path.read_text().splitlines())
    for header in lines:
        require(header == "ITEM: TIMESTEP", f"invalid frame in {path}")
        step = int(next(lines))
        require(next(lines) == "ITEM: NUMBER OF ATOMS", f"invalid count in {path}")
        require(int(next(lines)) == NATOMS, f"missing atoms in {path}")
        require(next(lines) == "ITEM: BOX BOUNDS pp pp pp", f"invalid box in {path}")
        for _ in range(3):
            bounds = [float(x) for x in next(lines).split()]
            require(bounds == [0, 4 * SPACING], f"changed box in {path}")
        require(
            next(lines) == "ITEM: ATOMS id xu yu zu fx fy fz",
            f"invalid columns in {path}",
        )
        atoms = []
        for expected in range(1, NATOMS + 1):
            row = next(lines).split()
            require(len(row) == 7 and int(row[0]) == expected, f"invalid ID in {path}")
            atoms.append([float(x) for x in row[1:]])
        require(
            all(math.isfinite(x) for row in atoms for x in row), f"nonfinite {path}"
        )
        frames.append((step, atoms))
    require(
        [step for step, _ in frames] == list(range(NSTEPS + 1)),
        f"missing frames in {path}",
    )
    return [atoms for _, atoms in frames]


def check_execution(directory, ranks, threads, biased, mode):
    log = (directory / "log.lammps").read_text()
    if mode == "mpi":
        grids = re.findall(r"(\d+) by (\d+) by (\d+) MPI processor grid", log)
        require(
            len(grids) == 1 and math.prod(map(int, grids[0])) == ranks,
            f"{directory.name}: LAMMPS did not use {ranks} MPI ranks",
        )
        require(
            re.search(rf"with {ranks} MPI tasks x {threads} OpenMP threads\b", log),
            f"{directory.name}: LAMMPS rank/thread count mismatch",
        )
    else:
        require(
            re.search(rf"using {threads} OpenMP thread\(s\) per MPI task\b", log),
            f"{directory.name}: serial LAMMPS OpenMP thread count mismatch",
        )
    require("lj/cut/omp" in log, f"{directory.name}: threaded pair style not active")
    if biased:
        log = (directory / "PLUMED.OUT").read_text()
        require(
            f"Running over {ranks} node" in log,
            f"{directory.name}: PLUMED rank count mismatch",
        )
        require(
            f"Number of threads: {threads}\n" in log,
            f"{directory.name}: PLUMED thread count mismatch",
        )


def invoke(args, ranks, threads, biased, directory):
    directory.mkdir()
    for name in ("argon.data", "in.lammps", "plumed.dat"):
        shutil.copy2(args.results / name, directory / name)
    command = []
    if args.mode == "mpi":
        command = [str(args.launcher)]
        if (args.launcher_kind or args.backend) == "openmpi":
            command += ["--oversubscribe", "--bind-to", "none"]
        command += ["-n", str(ranks)]
    command += [
        str(args.lmp),
        "-sf",
        "omp",
        "-pk",
        "omp",
        str(threads),
        "-var",
        "with_plumed",
        str(int(biased)),
        "-in",
        "in.lammps",
    ]
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    env = os.environ.copy()
    env.update(
        OMP_NUM_THREADS=str(threads),
        OMP_DYNAMIC="FALSE",
        PLUMED_NUM_THREADS=str(threads),
        PLUMED_KERNEL=str(args.kernel),
        MPICC="/usr/bin/false",
        MPICXX="/usr/bin/false",
    )
    if args.mode == "mpi":
        env.pop("PLUMED_MPI_WRAPPER_FORCE", None)
        env["PLUMED_MPI_WRAPPER"] = str(args.results / "must-not-load-provider")
    else:
        env.pop("PLUMED_MPI_WRAPPER", None)
        env.pop("PLUMED_MPI_WRAPPER_FORCE", None)
    with (directory / "stdout.log").open("w") as log:
        result = subprocess.run(
            command,
            cwd=directory,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=args.timeout,
        )
    if result.returncode:
        tail = "\n".join((directory / "stdout.log").read_text().splitlines()[-50:])
        raise RuntimeError(
            f"{directory.name}: LAMMPS failed ({result.returncode})\n{tail}"
        )


def check_results(directory):
    colvar = read_table(directory / "COLVAR")
    cvforces = read_table(directory / "CVFORCES")
    frames = read_dump(directory / "atoms.dump")
    require(
        len(colvar) == len(cvforces) == NSTEPS + 1,
        f"{directory.name}: missing PLUMED frames",
    )
    for step, (row, force, atoms) in enumerate(zip(colvar, cvforces, frames)):
        require(len(row) == 4 and len(force) == 2, "unexpected PLUMED columns")
        time, distance, bias, _energy = row
        require(
            abs(time - DT * step) < 1e-11 and abs(force[0] - time) < 1e-11,
            f"{directory.name}: incorrect time steps",
        )
        measured = math.dist(atoms[0][:3], atoms[1][:3])
        require(abs(distance - measured) < 2e-10, f"{directory.name}: incorrect CV")
        require(
            abs(bias - 0.5 * KAPPA * (distance - TARGET) ** 2) < 2e-10,
            f"{directory.name}: incorrect analytic bias",
        )
        require(
            abs(force[1] + KAPPA * (distance - TARGET)) < 2e-10,
            f"{directory.name}: incorrect generalized force",
        )
    # Native potential excludes the PLUMED bias (fix_modify energy no).
    log = (directory / "log.lammps").read_text()
    thermo = re.findall(
        r"^\s*(\d+)\s+([\deE.+-]+)\s+([\deE.+-]+)\s*$", log, re.MULTILINE
    )
    require(
        [int(row[0]) for row in thermo] == list(range(NSTEPS + 1)),
        f"{directory.name}: incomplete native potential energy",
    )
    energy_delta = max(
        abs(row[3] - float(native[2])) for row, native in zip(colvar, thermo)
    )
    require(
        energy_delta < 2e-9,
        f"{directory.name}: ENERGY differs from LAMMPS by {energy_delta}",
    )
    return colvar, frames


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("mpi", "serial"), default="mpi")
    parser.add_argument(
        "--backend", choices=("openmpi", "mpich", "serial"), required=True
    )
    parser.add_argument("--lmp", type=Path, required=True)
    parser.add_argument("--launcher", type=Path)
    parser.add_argument("--launcher-kind", choices=("openmpi", "mpich"))
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()
    if args.mode == "mpi" and args.launcher is None:
        parser.error("MPI mode requires --launcher")
    if (args.mode == "serial") != (args.backend == "serial"):
        parser.error("--mode serial and --backend serial must be used together")
    for name in ("lmp", "kernel"):
        path = getattr(args, name).resolve()
        require(path.is_file(), f"{name} not found: {path}")
        setattr(args, name, path)
    if args.launcher is not None:
        args.launcher = args.launcher.resolve()
        require(args.launcher.is_file(), f"launcher not found: {args.launcher}")
    args.results = args.results.resolve()
    if not args.check_only:
        args.results.mkdir(parents=True, exist_ok=False)
        for name in ("in.lammps", "plumed.dat"):
            shutil.copy2(ROOT / name, args.results / name)
        write_data(args.results / "argon.data")
    require(
        not (args.results / "must-not-load-provider").exists(),
        "provider sentinel exists",
    )
    print(f"Backend: {args.backend}; LAMMPS: {args.lmp}", flush=True)
    print(f"Results directory: {args.results}", flush=True)

    summary = {
        "backend": args.backend,
        "lmp": str(args.lmp),
        "kernel": str(args.kernel),
        "cases": {},
    }
    reference = None
    cases = MPI_CASES if args.mode == "mpi" else SERIAL_CASES
    for ranks, threads in cases:
        label = (
            f"r{ranks}-omp{threads}" if args.mode == "mpi" else f"serial-omp{threads}"
        )
        print(f"Checking {label}", flush=True)
        directory = args.results / label
        if not args.check_only:
            invoke(args, ranks, threads, True, directory)
        check_execution(directory, ranks, threads, True, args.mode)
        colvar, frames = check_results(directory)
        if reference is None:
            reference = colvar, frames
        delta_colvar = max(
            abs(a - b)
            for ref, row in zip(reference[0], colvar)
            for a, b in zip(ref, row)
        )
        delta_atoms = max(
            abs(a - b)
            for ref, frame in zip(reference[1], frames)
            for ref_atom, atom in zip(ref, frame)
            for a, b in zip(ref_atom, atom)
        )
        require(
            delta_colvar < 2e-9 and delta_atoms < 2e-8,
            f"{label}: rank/hybrid results differ",
        )
        print(
            f"{label}: max COLVAR/position/force deltas = {delta_colvar:.3g}/{delta_atoms:.3g}"
        )
        summary["cases"][label] = {
            "mode": args.mode,
            "threads": threads,
            "max_colvar_delta": delta_colvar,
            "max_atom_delta": delta_atoms,
        }
        if args.mode == "mpi":
            summary["cases"][label]["ranks"] = ranks

    baseline_ranks = 4 if args.mode == "mpi" else 1
    baseline_label = (
        f"baseline-r{baseline_ranks}-omp1"
        if args.mode == "mpi"
        else "baseline-serial-omp1"
    )
    baseline = args.results / baseline_label
    if not args.check_only:
        invoke(args, baseline_ranks, 1, False, baseline)
    check_execution(baseline, baseline_ranks, 1, False, args.mode)
    native = read_dump(baseline / "atoms.dump")[0]
    biased_label = f"r{baseline_ranks}-omp1" if args.mode == "mpi" else "serial-omp1"
    biased = read_dump(args.results / biased_label / "atoms.dump")[0]
    displacement = [b - a for a, b in zip(biased[0][:3], biased[1][:3])]
    distance = math.sqrt(sum(x * x for x in displacement))
    expected_force = [KAPPA * (distance - TARGET) * x / distance for x in displacement]
    force_error = 0.0
    for i, (without, with_bias) in enumerate(zip(native, biased)):
        sign = 1 if i == 0 else -1 if i == 1 else 0
        for axis in range(3):
            actual = with_bias[axis + 3] - without[axis + 3]
            force_error = max(force_error, abs(actual - sign * expected_force[axis]))
    require(force_error < 2e-9, f"incorrect atom force increment: {force_error}")
    summary["max_force_increment_error"] = force_error
    (args.results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if args.mode == "mpi":
        print(
            "PASS: injection, MPI/OpenMP counts, analytic bias/forces, native ENERGY, and rank/hybrid comparisons"
        )
    else:
        print(
            "PASS: serial execution, OpenMP counts, analytic bias/forces, native ENERGY, and thread comparisons"
        )


if __name__ == "__main__":
    main()
