#!/usr/bin/env python3
"""Run GROMACS 2026 PLUMED injection checks for tMPI, Open MPI, or MPICH."""

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
CASES = [(1, 1), (2, 1), (4, 1), (4, 2)]


def run(command, cwd, env, log_path, timeout):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        result = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    if result.returncode:
        tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-50:])
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n{tail}"
        )


def read_table(path):
    rows = [
        [float(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    if any(not math.isfinite(value) for row in rows for value in row):
        raise RuntimeError(f"non-finite value in {path}")
    return rows


def read_positions(path):
    lines = path.read_text().splitlines()
    coords = []
    for start in range(0, len(lines), 4):
        if start + 3 >= len(lines) or lines[start].strip() != "2":
            raise RuntimeError(f"malformed DUMPATOMS output: {path}")
        for line in lines[start + 2 : start + 4]:
            fields = line.split()
            coords.extend(float(value) for value in fields[1:4])
    if any(not math.isfinite(value) for value in coords):
        raise RuntimeError(f"non-finite coordinate in {path}")
    return coords


def first_frame_forces(gmx, traj, cwd, env, log_path, timeout):
    result = subprocess.run(
        [str(gmx), "dump", "-f", str(traj)],
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
        check=True,
    )
    log_path.write_text(result.stdout)
    try:
        first = result.stdout.split(" frame 0:", 1)[1].split(" frame 1:", 1)[0]
    except IndexError as exc:
        raise RuntimeError(f"no initial frame in {traj}") from exc
    forces = []
    for atom in (0, 1):
        match = re.search(rf"f\[\s*{atom}\]=\{{([^}}]+)\}}", first)
        if not match:
            raise RuntimeError(f"trajectory has no force for atom {atom + 1}: {traj}")
        forces.append([float(value) for value in match.group(1).split(",")])
    if any(len(force) != 3 for force in forces) or any(
        not math.isfinite(value) for force in forces for value in force
    ):
        raise RuntimeError(f"invalid force in {traj}")
    return forces


def invoke_mdrun(
    args, gmx, launcher, ranks, threads, directory, tpr, with_plumed, kernel
):
    directory.mkdir(parents=True, exist_ok=True)
    command = []
    if launcher is not None:
        command.append(str(launcher))
        if (args.launcher_kind or args.backend) == "openmpi":
            # CI can have fewer slots than ranks or hybrid worker threads.
            command.extend(["--oversubscribe", "--bind-to", "none"])
        command.extend(["-n", str(ranks)])
    command.extend([str(gmx), "mdrun", "-s", str(tpr)])
    if with_plumed:
        shutil.copy2(ROOT / "plumed.dat", directory / "plumed.dat")
        command.extend(["-plumed", "plumed.dat"])
    if args.backend == "tmpi":
        command.extend(["-ntmpi", str(ranks)])
    command.extend(["-ntomp", str(threads), "-pin", "off"])

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["PLUMED_NUM_THREADS"] = str(threads)
    env["PLUMED_KERNEL"] = str(kernel)
    env["MPICC"] = "/usr/bin/false"
    env["MPICXX"] = "/usr/bin/false"
    env.pop("PLUMED_MPI_WRAPPER_FORCE", None)
    if with_plumed:
        # Success proves GROMACS injection wins when ordinary environment
        # provider discovery points at a file that does not exist.
        env["PLUMED_MPI_WRAPPER"] = str(args.results / "must-not-load-provider")
    else:
        env.pop("PLUMED_MPI_WRAPPER", None)

    run(
        command,
        directory,
        env,
        directory / "run.log",
        args.timeout,
    )
    check_execution(args.backend, ranks, threads, directory, with_plumed)


def check_execution(backend, ranks, threads, directory, with_plumed=True):
    log = (directory / "run.log").read_text()
    rank_kind = "threads?" if backend == "tmpi" else "process(?:es)?"
    if not re.search(rf"Using {ranks} MPI {rank_kind}\b", log):
        raise RuntimeError(f"{directory.name}: GROMACS did not use {ranks} MPI ranks")
    if not re.search(rf"Using {threads} OpenMP threads?\b", log):
        raise RuntimeError(
            f"{directory.name}: GROMACS did not use {threads} OpenMP threads"
        )
    if with_plumed:
        log = (directory / "PLUMED.OUT").read_text()
        if f"Running over {ranks} node" not in log:
            raise RuntimeError(
                f"{directory.name}: PLUMED did not receive {ranks} MPI ranks"
            )
        if f"Number of threads: {threads}\n" not in log:
            raise RuntimeError(
                f"{directory.name}: PLUMED did not receive the thread count"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--backend", choices=("tmpi", "openmpi", "mpich"), required=True
    )
    parser.add_argument("--gmx", type=Path, required=True, help="GROMACS executable")
    parser.add_argument(
        "--gromacs-source",
        type=Path,
        help="pinned GROMACS source tree (required unless --check-only)",
    )
    parser.add_argument(
        "--launcher", type=Path, help="MPI launcher (required for real MPI)"
    )
    parser.add_argument(
        "--launcher-kind",
        choices=("openmpi", "mpich"),
        help="launcher implementation, if different from --backend",
    )
    parser.add_argument(
        "--kernel", type=Path, required=True, help="serial libplumedKernel"
    )
    parser.add_argument(
        "--results", type=Path, required=True, help="new output directory"
    )
    parser.add_argument(
        "--check-only", action="store_true", help="recheck existing results"
    )
    parser.add_argument("--timeout", type=int, default=180, help="seconds per mdrun")
    args = parser.parse_args()

    if not args.check_only:
        if args.gromacs_source is None:
            parser.error("--gromacs-source is required unless --check-only")
        argon = args.gromacs_source / "src/testutils/simulationdatabase/argon5832.gro"
        if not argon.is_file():
            parser.error(f"GROMACS argon fixture not found: {argon}")

    gmx = args.gmx.resolve()
    if not gmx.is_file():
        raise SystemExit(f"GROMACS executable not found: {gmx}")

    launcher = None
    if args.backend != "tmpi":
        if args.launcher is None:
            parser.error("--launcher is required for Open MPI and MPICH")
        launcher = args.launcher.resolve()
        if not launcher.is_file():
            raise SystemExit(f"MPI launcher not found: {launcher}")
    elif args.launcher is not None or args.launcher_kind is not None:
        parser.error("thread-MPI must not use an external MPI launcher")

    kernel = args.kernel.resolve()
    if not kernel.is_file():
        raise SystemExit(f"serial PLUMED kernel not found: {kernel}")

    # grompp is a serial preprocessing step; run it once with this backend's
    # executable and keep the generated TPR in the dedicated test workspace.
    prep = args.results.resolve()
    args.results = prep
    if not args.check_only:
        prep.mkdir(parents=True, exist_ok=False)
    if (prep / "must-not-load-provider").exists():
        raise SystemExit("injection sentinel must not exist")
    print(f"Backend: {args.backend}; GROMACS executable: {gmx}", flush=True)
    print(f"Results directory: {prep}", flush=True)
    tpr = prep / "base.tpr"
    base_env = os.environ.copy()
    if not args.check_only:
        # Reuse GROMACS's own fixture without duplicating its large data file.
        shutil.copy2(argon, prep / "argon5832.gro")
        for name in ("argon5832.top", "run.mdp", "plumed.dat"):
            shutil.copy2(ROOT / name, prep / name)
        run(
            [
                str(gmx),
                "grompp",
                "-f",
                "run.mdp",
                "-c",
                "argon5832.gro",
                "-p",
                "argon5832.top",
                "-o",
                str(tpr),
            ],
            prep,
            base_env,
            prep / "grompp.log",
            args.timeout,
        )
    # mdrun reads the TPR from this backend-specific test directory.
    outputs = {}
    for ranks, threads in CASES:
        label = f"r{ranks}-omp{threads}"
        print(f"Running {label} with PLUMED", flush=True)
        with_dir = prep / f"plumed-{label}"
        outputs[label] = with_dir
        if not args.check_only:
            invoke_mdrun(
                args, gmx, launcher, ranks, threads, with_dir, tpr, True, kernel
            )
        else:
            check_execution(args.backend, ranks, threads, with_dir)
        if label == "r4-omp1" and not args.check_only:
            baseline = prep / "baseline-r4-omp1"
            print("Running four-rank no-PLUMED force baseline", flush=True)
            invoke_mdrun(
                args, gmx, launcher, ranks, threads, baseline, tpr, False, kernel
            )

    summary = {
        "backend": args.backend,
        "gmx": str(gmx),
        "kernel": str(kernel),
        "cases": {},
    }
    reference = None
    reference_xyz = None
    for label, directory in outputs.items():
        colvar = read_table(directory / "COLVAR")
        cvforces = read_table(directory / "CVFORCES")
        xyz = read_positions(directory / "atoms.xyz")
        if any(len(row) != 4 for row in colvar) or any(
            len(row) != 2 for row in cvforces
        ):
            raise RuntimeError(f"{label}: unexpected output columns")
        if len(colvar) != 21 or len(cvforces) != 21 or len(xyz) != 126:
            raise RuntimeError(f"{label}: expected 21 output frames")
        for step, (row, force_row) in enumerate(zip(colvar, cvforces)):
            time, distance, bias, _energy = row
            force_time, cv_force = force_row
            expected_bias = 50.0 * (distance - 0.3) ** 2
            expected_cv_force = -100.0 * (distance - 0.3)
            if abs(time - force_time) > 1e-12:
                raise RuntimeError(f"{label}: COLVAR/CVFORCES time mismatch")
            if abs(time - 0.0005 * step) > 1e-10:
                raise RuntimeError(f"{label}: missing or repeated time step")
            if abs(bias - expected_bias) > 2e-8:
                raise RuntimeError(f"{label}: restraint bias is not analytic")
            if abs(cv_force - expected_cv_force) > 2e-8:
                raise RuntimeError(f"{label}: generalized force is not analytic")
        if reference is None:
            reference, reference_xyz = colvar, xyz
            summary["cases"][label] = {"reference": True}
        else:
            delta_distance = max(abs(a[1] - b[1]) for a, b in zip(reference, colvar))
            delta_bias = max(abs(a[2] - b[2]) for a, b in zip(reference, colvar))
            delta_energy = max(abs(a[3] - b[3]) for a, b in zip(reference, colvar))
            delta_xyz = max(abs(a - b) for a, b in zip(reference_xyz, xyz))
            if delta_distance > 2e-6 or delta_bias > 2e-4 or delta_xyz > 2e-5:
                raise RuntimeError(f"{label}: positions/CV differ from 1-rank run")
            if delta_energy > 5e-6 * max(1.0, abs(reference[0][3])):
                raise RuntimeError(f"{label}: parallel energy reduction diverged")
            print(
                f"{label}: max distance/bias/energy/position deltas = "
                f"{delta_distance:.3g}/{delta_bias:.3g}/{delta_energy:.3g}/{delta_xyz:.3g}"
            )
            summary["cases"][label] = {
                "max_distance_delta": delta_distance,
                "max_bias_delta": delta_bias,
                "max_energy_delta": delta_energy,
                "max_position_delta": delta_xyz,
            }

    biased_forces = first_frame_forces(
        gmx,
        outputs["r4-omp1"] / "traj.trr",
        prep,
        os.environ.copy(),
        prep / "biased-first-frame.dump",
        args.timeout,
    )
    baseline_forces = first_frame_forces(
        gmx,
        prep / "baseline-r4-omp1/traj.trr",
        prep,
        os.environ.copy(),
        prep / "baseline-first-frame.dump",
        args.timeout,
    )
    xyz_rows = read_positions(outputs["r4-omp1"] / "atoms.xyz")
    delta = [xyz_rows[3 + i] - xyz_rows[i] for i in range(3)]
    distance = math.sqrt(sum(value * value for value in delta))
    cv_force = -100.0 * (distance - 0.3)
    expected_atom_1 = [-cv_force * value / distance for value in delta]
    expected_atom_2 = [-value for value in expected_atom_1]
    for atom, expected in enumerate((expected_atom_1, expected_atom_2)):
        actual = [biased_forces[atom][i] - baseline_forces[atom][i] for i in range(3)]
        if any(abs(a - b) > 0.02 for a, b in zip(actual, expected)):
            raise RuntimeError(
                f"PLUMED force was not added to GROMACS atom {atom + 1}: "
                f"observed {actual}, expected {expected}"
            )
    print("PASS: analytic bias, CV forces, atom forces, and rank/hybrid comparisons")
    (prep / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Results directory: {prep}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, subprocess.SubprocessError, RuntimeError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error
