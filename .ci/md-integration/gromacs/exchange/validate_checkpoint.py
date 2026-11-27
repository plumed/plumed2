"""Check common-step replica checkpoints with and without shared MAXENT.

Four two-particle replicas accept swap 0/1 and reject swap 2/3 at step 200.
With replica-local checkpoint eligibility, the pairs checkpoint at 201/250;
MAXENT's broadcast at 250 then deadlocks against the native checkpoint barrier.
No transfer/group extensions or third-party Python packages are needed.
"""

import argparse
import json
import math
import os
import re
import shlex
import shutil
import signal
import subprocess
import time
from pathlib import Path

TOP = """[ defaults ]
1 2 yes 0.5 0.8333333333
[ atomtypes ]
X 1000.0 0.0 A 0.3 0.0
[ moleculetype ]
DUM 0
[ atoms ]
1 X 1 DUM X1 1 0.0 1000.0
2 X 1 DUM X2 1 0.0 1000.0
[ system ]
Checkpoint collective test
[ molecules ]
DUM 1
"""
MDP = """integrator = md
dt = .001
nsteps = 400
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = 0
nstenergy = 0
nstlog = 0
nstcalcenergy = 10
cutoff-scheme = Verlet
nstlist = 50
verlet-buffer-tolerance = -1
rlist = 1.1
rcoulomb = 1.0
rvdw = 1.0
coulombtype = Cut-off
vdwtype = Cut-off
constraints = none
pbc = xyz
tcoupl = V-rescale
tc-grps = System
tau-t = .1
ref-t = 300
pcoupl = no
comm-mode = Linear
nstcomm = 10
gen-vel = yes
gen-temp = 300
"""


def run(command, directory, env, log, timeout):
    """Bound both initial and restart runs, and reap the entire MPI process group."""
    (directory / (log + ".command.json")).write_text(json.dumps(command) + "\n")
    started = time.monotonic()
    timed_out = False
    with (directory / log).open("w") as output:
        child = subprocess.Popen(
            command,
            cwd=directory,
            env=env,
            stdout=output,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            code = child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(child.pid, signal.SIGTERM)
            try:
                code = child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                code = child.wait()
    return {
        "exit_code": code,
        "timed_out": timed_out,
        "wall_seconds": time.monotonic() - started,
    }


def rows(path):
    return [
        list(map(float, line.split()))
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def prepare(args, directory, shared, env):
    for replica in range(4):
        rep = directory / f"R{replica}"
        rep.mkdir()
        (rep / "topol.top").write_text(TOP)
        (rep / "test.mdp").write_text(
            MDP + f"gen-seed = {42 + replica}\nld-seed = {754 + replica}\n"
        )
        distance = 0.8 if replica == 3 else 0.3
        coordinates = "Checkpoint test\n    2\n"
        for index, name, x in [(1, "X1", 1.0), (2, "X2", 1.0 + distance)]:
            coordinates += (
                f"{1:5d}{'DUM':<5}{name:>5}{index:5d}{x:8.3f}{1.0:8.3f}{1.0:8.3f}\n"
            )
        (rep / "conf.gro").write_text(coordinates + "4.0 4.0 4.0\n")
        result = run(
            [
                args.gmx,
                "grompp",
                "-f",
                "test.mdp",
                "-c",
                "conf.gro",
                "-p",
                "topol.top",
                "-o",
                "md.tpr",
            ],
            rep,
            env,
            "grompp.log",
            args.timeout,
        )
        if result["timed_out"] or result["exit_code"] != 0:
            raise RuntimeError(f"grompp failed: {rep / 'grompp.log'}: {result}")
        lines = [
            "CHECKPOINT FILE=plumed.cpt",
            "d: DISTANCE ATOMS=1,2 NOPBC",
            f"b: RESTRAINT ARG=d AT={distance} KAPPA={10000 if replica >= 2 else 0}",
            "PRINT ARG=d FILE=trace STRIDE=1 FMT=%.17g",
            "FLUSH STRIDE=1",
        ]
        if shared:
            # Zero force weights isolate the collective from any bias-force changes.
            lines += [
                (
                    "m: MAXENT ARG=d TYPE=EQUAL AT=.3 KAPPA=.001 TAU=10 "
                    "PACE=50 PRINT_STRIDE=50 FMT=%.17g TEMP=300 "
                    "LEARN_REPLICA=0 APPLY_WEIGHTS=0,0,0,0"
                )
            ]
        (rep / "plumed.dat").write_text("\n".join(lines) + "\n")


def check_case(args, root, mode, shared, env):
    ranks = 2 if mode == "dd" else 1
    name = mode + ("-maxent" if shared else "-control")
    directory = root / name
    directory.mkdir()
    result = {"name": name, "ranks_per_replica": ranks, "shared_maxent": shared}
    try:
        prepare(args, directory, shared, env)
        command = [part.format(n=4 * ranks) for part in shlex.split(args.launcher)]
        command += [
            args.gmx,
            "mdrun",
            "-multidir",
            "R0",
            "R1",
            "R2",
            "R3",
            "-deffnm",
            "md",
            "-plumed",
            "plumed.dat",
            "-replex",
            "200",
            "-reseed",
            "20261123",
            "-cpt",
            "0",
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
            "-notunepme",
        ]
        if mode == "hrex":
            command += ["-hrex"]
        result.update(run(command, directory, env, "mdrun.log", args.timeout))
        checkpoints, last_trace, exchanges = {}, {}, {}
        for replica in range(4):
            rep = directory / f"R{replica}"
            log = (rep / "md.log").read_text()
            rank_count = re.search(r"Using (\d+) MPI process", log)
            if rank_count is None or int(rank_count[1]) != ranks:
                raise AssertionError(f"Expected {ranks} MPI ranks per replica in {rep}")
            if mode == "dd" and "Domain decomposition grid" not in log:
                raise AssertionError(f"Domain decomposition was not exercised in {rep}")
            checkpoints[str(replica)] = [
                int(s) for s in re.findall(r"Writing checkpoint, step (\d+)", log)
            ]
            trace = rows(rep / f"trace.{replica}")
            last_trace[str(replica)] = trace[-1][0] if trace else None
            exchanges[str(replica)] = [
                line.split() for line in log.splitlines() if line.startswith("Repl ex")
            ]
        result.update(
            checkpoint_steps=checkpoints,
            last_trace_time_ps=last_trace,
            exchanges=exchanges,
        )
        if result["timed_out"] or result["exit_code"] != 0:
            raise RuntimeError(f"MD failed or hung: {result}")
        if any(steps != [200, 400] for steps in checkpoints.values()):
            raise AssertionError(
                f"Expected common checkpoint steps 200,400: {checkpoints}"
            )
        expected_swap = [["Repl", "ex", "0", "x", "1", "2", "3"]]
        if any(value != expected_swap for value in exchanges.values()):
            raise AssertionError(f"Fixture must accept only swap 0/1: {exchanges}")
        names = ["trace"] + (["m.LAGMULT"] if shared else [])
        reference = {}
        for replica in range(4):
            rep = directory / f"R{replica}"
            trace = rows(rep / f"trace.{replica}")
            if len(trace) != 401 or any(
                not math.isclose(row[0], step * 0.001, abs_tol=1e-14)
                or not all(math.isfinite(value) for value in row)
                for step, row in enumerate(trace)
            ):
                raise AssertionError(f"Incomplete or invalid trace in {rep}")
            if shared:
                multipliers = rows(rep / f"m.LAGMULT.{replica}")
                if (
                    len(multipliers) != 9
                    or not all(math.isfinite(v) for row in multipliers for v in row)
                    or not any(row[1] != 0 for row in multipliers)
                ):
                    raise AssertionError(f"Missing MAXENT learning in {rep}")
            for name in names:
                reference[replica, name] = (rep / f"{name}.{replica}").read_bytes()
            # The previous set is the INTERMEDIATE step-200 checkpoint, not the final one.
            shutil.copyfile(rep / "md_prev.cpt", rep / "md.cpt")
            shutil.copyfile(
                rep / f"plumed_prev.{replica}.cpt", rep / f"plumed.{replica}.cpt"
            )
        if shared and any(
            reference[replica, "m.LAGMULT"] != reference[0, "m.LAGMULT"]
            for replica in range(1, 4)
        ):
            raise AssertionError("MAXENT multipliers differ between replicas")
        restart = run(
            command + ["-cpi", "md.cpt", "-append"],
            directory,
            env,
            "restart.log",
            args.timeout,
        )
        result["restart"] = restart
        if restart["timed_out"] or restart["exit_code"] != 0:
            raise RuntimeError(f"Intermediate restart failed: {restart}")
        identical = all(
            (directory / f"R{replica}/{name}.{replica}").read_bytes() == data
            for (replica, name), data in reference.items()
        )
        result["restart_trace_and_coupling_bitwise_identical"] = identical
        if not identical:
            raise AssertionError(
                "Intermediate restart changed trace or MAXENT multipliers"
            )
        result["status"] = "pass"
    except (AssertionError, OSError, RuntimeError, ValueError) as error:
        result.update(status="fail", error=str(error))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--launcher", required=True, help="MPI command containing {n}")
    parser.add_argument("--output", required=True, help="Fresh results directory")
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--case", action="append", choices=["ordinary", "dd", "hrex"])
    args = parser.parse_args()
    if "{n}" not in args.launcher:
        parser.error("--launcher must contain {n}")
    args.gmx = str(Path(args.gmx).resolve())
    args.kernel = str(Path(args.kernel).resolve())
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / "arguments.json").write_text(json.dumps(vars(args), indent=2) + "\n")
    env = os.environ.copy()
    env.update(
        PLUMED_KERNEL=args.kernel,
        OMP_NUM_THREADS="1",
        PLUMED_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        GMX_MAXBACKUP="-1",
        GMX_DISABLE_MODULAR_SIMULATOR="1",
    )
    for variable in ("PLUMED_MPI_WRAPPER", "PLUMED_MPI_WRAPPER_FORCE"):
        env.pop(variable, None)
    results = []
    for mode in args.case or ["ordinary", "dd", "hrex"]:
        for shared in (False, True):
            result = check_case(args, root, mode, shared, env)
            results.append(result)
            (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
            print(result["name"], result["status"], result.get("error", ""), flush=True)
    return int(any(result["status"] != "pass" for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
