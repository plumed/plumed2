"""HREX automatic-list synchronization with independent LJ/force/virial oracles.

Use a tight valid buffer tolerance to expose automatic-radius differences in a
64-atom system. All physical cutoffs stay at 1 nm. No project-specific groups,
transfer actions or hard-coded replacement buffer radius are used.
"""

import argparse
import json
import math
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path

from validate import R, table
from validate_checkpoint import run

DT = 0.001
STEPS = 40
NAMES = [
    "auto-nvt",
    "auto-npt",
    "auto-dd",
    "homogeneous",
    "manual",
    "no-hrex",
    "cutoff-mismatch",
    "box-incompatible",
]


def invoke(args, command, directory, env, log, stdin=None):
    if stdin is None and command[1] != "dump":
        result = run(command, directory, env, log, args.timeout)
        if result["timed_out"] or result["exit_code"]:
            raise RuntimeError(f"{directory / log}: {result}")
    else:
        with (
            (directory / log).open("w") as output,
            (directory / (log + ".stderr")).open("w") as error,
        ):
            subprocess.run(
                command,
                cwd=directory,
                env=env,
                input=stdin,
                text=True,
                stdout=output,
                stderr=error,
                timeout=args.timeout,
                check=True,
            )


def setup(args, directory, case, env):
    boxes = [2.4, 3.2, 2.4, 3.2] if case == "box-incompatible" else [3.2] * 4
    if case == "no-hrex":
        boxes = [3.2] * 4
    epsilons = [0.1] * 4 if case in ("homogeneous", "no-hrex") else [0.2, 0.1, 0.2, 0.1]
    automatic = case.startswith("auto") or case == "homogeneous"
    targets = []
    for replica, (box, epsilon) in enumerate(zip(boxes, epsilons)):
        rep = directory / f"rep{replica}"
        rep.mkdir()
        top = f"""[ defaults ]
1 2 no
[ atomtypes ]
Ar 39.948 0.0 A 0.34 {epsilon}
[ moleculetype ]
AR 1
[ atoms ]
1 Ar 1 AR AR 1 0.0 39.948
[ system ]
HREX buffer test
[ molecules ]
AR 64
"""
        (rep / "topol.top").write_text(top)
        coords = [
            [box * (x + 0.5) / 4, box * (y + 0.5) / 4, box * (z + 0.5) / 4]
            for x in range(4)
            for y in range(4)
            for z in range(4)
        ]
        # Distinct configurations force accepted/rejected trials without changing
        # the NVT box; a displaced atom also avoids a symmetry-zero force oracle.
        coords[0][0] += 0.03 if replica % 2 == 0 or case == "no-hrex" else 0.23
        coords = [[round(value, 3) for value in row] for row in coords]
        targets.append(math.dist(coords[0], coords[1]))
        gro = "Argon grid\n   64\n"
        for atom, xyz in enumerate(coords, 1):
            gro += (
                f"{atom:5d}{'AR':<5}{'AR':>5}{atom:5d}"
                + "".join(f"{x:8.3f}" for x in xyz)
                + "\n"
            )
        (rep / "conf.gro").write_text(gro + f"{box} {box} {box}\n")
        cutoff = 0.9 if case == "cutoff-mismatch" and replica == 3 else 1.0
        radius = (
            [1.05, 1.35, 1.05, 1.35][replica]
            if case == "box-incompatible"
            else [1.1, 1.2, 1.15, 1.25][replica]
        )
        mdp = f"""integrator = md
nsteps = {STEPS}
dt = {DT}
nstlist = 20
cutoff-scheme = Verlet
rlist = {radius}
verlet-buffer-tolerance = {1e-6 if automatic else -1}
coulombtype = Cut-off
rcoulomb = {cutoff}
vdwtype = Cut-off
vdw-modifier = None
rvdw = {cutoff}
pbc = xyz
tcoupl = v-rescale
tc-grps = System
tau-t = .1
ref-t = 300
pcoupl = {"C-rescale" if case == "auto-npt" else "no"}
constraints = none
gen-vel = yes
gen-temp = 300
gen-seed = {12345 + replica}
ld-seed = {23456 + replica}
nstcomm = 10
nstenergy = 1
nstlog = 0
nstcalcenergy = 1
nstxout = 1
nstvout = 0
nstfout = 1
"""
        if case == "auto-npt":
            mdp += "pcoupltype = isotropic\ntau-p = 5\nref-p = 1\ncompressibility = 4.5e-5\nnstpcouple = 10\n"
        (rep / "md.mdp").write_text(mdp)
        kappa = 0 if replica < 2 else 100000
        (rep / "plumed.dat").write_text(
            "CHECKPOINT FILE=plumed.cpt\nd: DISTANCE ATOMS=1,2\ne: ENERGY\n"
            f"b: RESTRAINT ARG=d AT={targets[-1]:.17g} KAPPA={kappa}\n"
            "PRINT ARG=d,e,b.bias FILE=COLVAR STRIDE=1 FMT=%.17g\nFLUSH STRIDE=1\n"
        )
        invoke(
            args,
            [
                args.gmx,
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
            rep,
            env,
            "grompp.log",
        )
    return epsilons, targets


def frames(text):
    result = {}
    for block in re.split(r"^.* frame \d+:\s*$", text, flags=re.M)[1:]:
        step = int(re.search(r"\bstep=\s*(\d+)", block)[1])
        values = {}
        for name in ("box", "x", "f"):
            values[name] = [
                [float(x) for x in m.split(",")]
                for m in re.findall(rf"\b{name}\[\s*\d+\]\s*=\{{([^}}]+)\}}", block)
            ]
        if len(values["x"]) != 64 or len(values["f"]) != 64 or len(values["box"]) != 3:
            raise AssertionError(f"Incomplete trajectory frame {step}")
        result[step] = values
    if sorted(result) != list(range(STEPS + 1)):
        raise AssertionError("Missing trajectory steps")
    return result


def lj(frame, epsilon):
    """Direct O(N^2) minimum-image LJ sum; no engine-reported energy/force inputs."""
    xyz = frame["x"]
    box = [frame["box"][i][i] for i in range(3)]
    if any(
        abs(frame["box"][i][j]) > 1e-10 for i in range(3) for j in range(3) if i != j
    ):
        raise AssertionError("Oracle expects orthorhombic boxes")
    force = [[0.0] * 3 for _ in xyz]
    virial = [[0.0] * 3 for _ in range(3)]
    energy = 0.0
    for i in range(len(xyz)):
        for j in range(i):
            delta = [xyz[i][k] - xyz[j][k] for k in range(3)]
            delta = [d - round(d / size) * size for d, size in zip(delta, box)]
            r2 = sum(d * d for d in delta)
            if r2 >= 1.0:
                continue
            z = (0.34 * 0.34 / r2) ** 3
            energy += 4 * epsilon * (z * z - z)
            coefficient = 24 * epsilon * (2 * z * z - z) / r2
            f = [coefficient * d for d in delta]
            for a in range(3):
                force[i][a] += f[a]
                force[j][a] -= f[a]
                for b in range(3):
                    virial[a][b] -= 0.5 * delta[a] * f[b]
    return energy, force, virial


def verify(args, directory, epsilons, targets, env):
    trajectories, rows = [], []
    maxima = dict(energy=0.0, force=0.0, virial=0.0, bias=0.0)
    for replica in range(4):
        rep = directory / f"rep{replica}"
        colvar = table(rep / f"COLVAR.{replica}")
        if len(colvar) != STEPS + 1:
            raise AssertionError("Missing/duplicate ordinary PLUMED steps")
        rows.append(colvar)
        invoke(args, [args.gmx, "dump", "-f", "run.trr"], rep, env, "trajectory.dump")
        trajectory = frames((rep / "trajectory.dump").read_text())
        trajectories.append(trajectory)
        components = ["Vir-" + a + b for a in "XYZ" for b in "XYZ"]
        invoke(
            args,
            [args.gmx, "energy", "-f", "run.edr", "-o", "virial.xvg", "-dp"],
            rep,
            env,
            "energy.log",
            stdin="\n".join(components + ["0", ""]),
        )
        native_virials = [
            [float(x) for x in line.split()]
            for line in (rep / "virial.xvg").read_text().splitlines()
            if line and line[0] not in "#@"
        ]
        if len(native_virials) != STEPS + 1 or any(
            len(r) != 10 for r in native_virials
        ):
            raise AssertionError("Missing native virial tensor")
        for step, row in enumerate(colvar):
            if not math.isclose(row["time"], step * DT, abs_tol=1e-12):
                raise AssertionError("Invalid PLUMED times")
            frame = trajectory[step]
            energy, force, virial = lj(frame, epsilons[replica])
            kappa = 0 if replica < 2 else 100000
            bias = 0.5 * kappa * (row["d"] - targets[replica]) ** 2
            delta = [frame["x"][0][k] - frame["x"][1][k] for k in range(3)]
            box = [frame["box"][k][k] for k in range(3)]
            delta = [d - round(d / size) * size for d, size in zip(delta, box)]
            distance = math.sqrt(sum(x * x for x in delta))
            if abs(distance - row["d"]) > 3e-6:
                raise AssertionError("Distance differs from native coordinates")
            bf = [-kappa * (row["d"] - targets[replica]) * d / distance for d in delta]
            for a in range(3):
                force[0][a] += bf[a]
                force[1][a] -= bf[a]
                for b in range(3):
                    virial[a][b] -= 0.5 * delta[a] * bf[b]
            maxima["energy"] = max(maxima["energy"], abs(row["e"] - energy))
            maxima["bias"] = max(maxima["bias"], abs(row["b.bias"] - bias))
            maxima["force"] = max(
                maxima["force"],
                max(
                    abs(a - b) for f, g in zip(force, frame["f"]) for a, b in zip(f, g)
                ),
            )
            maxima["virial"] = max(
                maxima["virial"],
                max(
                    abs(virial[a][b] - native_virials[step][1 + 3 * a + b])
                    for a in range(3)
                    for b in range(3)
                ),
            )
    budgets = dict(energy=5e-4, bias=2e-8, force=1e-2, virial=1e-2)
    for field, maximum in maxima.items():
        if maximum > budgets[field] or not math.isfinite(maximum):
            raise AssertionError(
                f"{field} oracle error {maximum} exceeds {budgets[field]}"
            )
    events = {}
    for replica in range(4):
        log = (directory / f"rep{replica}/run.log").read_text()
        for m in re.finditer(
            r"PLUMED_GREX step=(\d+) a=(\d+) b=(\d+) correction=(\S+) delta=(\S+)", log
        ):
            events[tuple(map(int, m.group(1, 2, 3)))] = tuple(map(float, m.group(4, 5)))
    if sorted(events) != [(20, 0, 1), (20, 2, 3)]:
        raise AssertionError(f"Unexpected exchange events {events}")
    exchanges = []
    for (step, a, b), (correction, total) in sorted(events.items()):
        ua = lj(trajectories[a][step], 1.0)[0]
        ub = lj(trajectories[b][step], 1.0)[0]
        physical = (epsilons[a] - epsilons[b]) * (ub - ua)
        bias = 0.0
        for receiver, donor in ((a, b), (b, a)):
            kappa = 0 if receiver < 2 else 100000
            bias += (
                0.5
                * kappa
                * (
                    (rows[donor][step]["d"] - targets[receiver]) ** 2
                    - (rows[receiver][step]["d"] - targets[receiver]) ** 2
                )
            )
        expected = (physical + bias) / (R * 300)
        if max(abs(correction - expected), abs(total - expected)) > 1e-3 + 2e-5 * abs(
            expected
        ):
            raise AssertionError(
                f"Cross-energy mismatch {a}/{b}: {correction}, {total} vs {expected}"
            )
        log = (directory / f"rep{a}/run.log").read_text()
        segment = log.split("Replica exchange at step 20 ")[-1]
        accepted = bool(re.search(rf"(?<!\d){a}\s+x\s+{b}(?!\d)", segment))
        if directory.name != "no-hrex" and accepted != (a == 0):
            raise AssertionError(
                f"Expected acceptance only for 0/1: pair {a}/{b} {accepted}"
            )
        for receiver, donor in ((a, b if accepted else a), (b, a if accepted else b)):
            if abs(rows[receiver][21]["d"] - rows[donor][20]["d"]) > 0.01:
                raise AssertionError("Wrong configuration moved after trial")
        exchanges.append(
            dict(a=a, b=b, expected_delta=expected, delta=total, accepted=accepted)
        )
    return dict(max_errors=maxima, budgets=budgets, exchanges=exchanges), rows


def check_lists(directory, hrex, phases):
    logs = [(directory / f"rep{i}/run.log").read_text() for i in range(4)]
    matches = [
        [
            tuple(map(float, m))
            for m in re.findall(
                r"HREX neighbor-list cutoff: local (\S+) nm, common (\S+) nm", log
            )
        ]
        for log in logs
    ]
    if not hrex:
        if any(matches):
            raise AssertionError(
                "Ordinary exchange unexpectedly synchronized HREX lists"
            )
        return []
    if any(len(rows) != phases for rows in matches):
        raise AssertionError(
            f"Missing startup/restart common-list diagnostics: {matches}"
        )
    for phase in range(phases):
        local = [rows[phase][0] for rows in matches]
        common = [rows[phase][1] for rows in matches]
        if len(set(common)) != 1 or common[0] != max(local):
            raise AssertionError(
                f"Invalid common conservative radius: {local} -> {common}"
            )
    if any(
        "Multi-checking HREX neighbor-list cutoff ... OK" not in log for log in logs
    ):
        raise AssertionError("Existing HREX equality guard was not exercised")
    return matches


def case_run(args, root, case, env):
    directory = root / case
    directory.mkdir()
    result = dict(name=case)
    try:
        epsilons, targets = setup(args, directory, case, env)
        ranks = 2 if case == "auto-dd" else 1
        hrex = case != "no-hrex"
        command = [part.format(n=4 * ranks) for part in shlex.split(args.launcher)] + [
            args.gmx,
            "mdrun",
            "-multidir",
            "rep0",
            "rep1",
            "rep2",
            "rep3",
            "-deffnm",
            "run",
            "-s",
            "topol.tpr",
            "-plumed",
            "plumed.dat",
            "-replex",
            "20",
            "-reseed",
            "20261123",
            "-notunepme",
            "-ntomp",
            "1",
            "-pin",
            "off",
            "-nb",
            "cpu",
            "-pme",
            "cpu",
            "-update",
            "cpu",
            "-cpt",
            "0",
        ]
        if hrex:
            command += ["-hrex"]
        result["run"] = run(command, directory, env, "mdrun.log", args.timeout)
        logs = (
            (directory / "mdrun.log").read_text()
            + "\n"
            + "\n".join(p.read_text() for p in directory.glob("rep*/run.log"))
        )
        if case in ("cutoff-mismatch", "box-incompatible"):
            expected = (
                "HREX Coulomb cutoff is not equal"
                if case == "cutoff-mismatch"
                else "common HREX neighbor-list cutoff"
            )
            if (
                result["run"]["timed_out"]
                or result["run"]["exit_code"] == 0
                or expected not in logs
            ):
                raise AssertionError(f"Expected rejection missing: {expected}")
            result["expected_rejection"] = expected
        else:
            if result["run"]["timed_out"] or result["run"]["exit_code"]:
                raise RuntimeError(
                    f"MD refused/hung; inspect {directory / 'mdrun.log'}"
                )
            for replica in range(4):
                log = (directory / f"rep{replica}/run.log").read_text()
                if re.findall(r"Writing checkpoint, step (\d+)", log) != ["20", "40"]:
                    raise AssertionError("Checkpoint synchronization changed")
                if f"Using {ranks} MPI process" not in log:
                    raise AssertionError("Wrong MPI rank count")
                if ranks > 1 and "Domain decomposition grid" not in log:
                    raise AssertionError("DD not exercised")
            lists = check_lists(directory, hrex, 1)
            if case == "auto-nvt" and len({row[0][0] for row in lists}) < 2:
                raise AssertionError(
                    "Auto-buffer fixture did not select distinct local radii"
                )
            result["oracles"], reference = verify(
                args, directory, epsilons, targets, env
            )
            for replica in range(4):
                rep = directory / f"rep{replica}"
                shutil.copyfile(rep / "run_prev.cpt", rep / "run.cpt")
                shutil.copyfile(
                    rep / f"plumed_prev.{replica}.cpt", rep / f"plumed.{replica}.cpt"
                )
            invoke(
                args,
                command + ["-cpi", "run.cpt", "-append"],
                directory,
                env,
                "restart.log",
            )
            result["radii"] = check_lists(directory, hrex, 2)
            result["restart_oracles"], resumed = verify(
                args, directory, epsilons, targets, env
            )
            maxima = {
                field: max(
                    abs(a[field] - b[field])
                    for x, y in zip(reference, resumed)
                    for a, b in zip(x, y)
                )
                for field in ("d", "e", "b.bias")
            }
            for field, tolerance in (("d", 2e-5), ("e", 5e-4), ("b.bias", 2e-2)):
                if maxima[field] > tolerance:
                    raise AssertionError(f"Restart changed {field} by {maxima[field]}")
            result["restart_max_errors"] = maxima
        result["status"] = "pass"
    except (
        AssertionError,
        RuntimeError,
        OSError,
        ValueError,
        KeyError,
        IndexError,
        subprocess.SubprocessError,
    ) as error:
        result.update(status="fail", error=str(error))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--launcher", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--case", choices=NAMES, action="append")
    parser.add_argument("--timeout", type=float, default=120)
    args = parser.parse_args()
    args.gmx = str(Path(args.gmx).resolve())
    args.kernel = str(Path(args.kernel).resolve())
    if "{n}" not in args.launcher:
        parser.error("launcher must contain {n}")
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(
        os.environ,
        PLUMED_KERNEL=args.kernel,
        OMP_NUM_THREADS="1",
        PLUMED_NUM_THREADS="1",
        GMX_MAXBACKUP="-1",
        GMX_PRINT_LONGFORMAT="1",
        GMX_DISABLE_MODULAR_SIMULATOR="1",
    )
    for key in ("PLUMED_MPI_WRAPPER", "PLUMED_MPI_WRAPPER_FORCE"):
        env.pop(key, None)
    (root / "arguments.json").write_text(json.dumps(vars(args), indent=2) + "\n")
    results = []
    for case in args.case or NAMES:
        result = case_run(args, root, case, env)
        results.append(result)
        (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
        print(case, result["status"], result.get("error", ""), flush=True)
    return int(any(r["status"] != "pass" for r in results))


if __name__ == "__main__":
    raise SystemExit(main())
