"""Validate alchemical metadynamics trials and lambda-state sampling independently."""

import argparse
import importlib.util
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location(
    "lambda_suite", Path(__file__).with_name("validate_lambda.py")
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def cases_for(quick=False, skip_population=False):
    cases = []
    for kind in (
        "lambda",
        "energy",
        "lambda-energy",
        "metad",
        "metad2d",
        "metad-energy",
    ):
        for ranks in [1, 2]:
            for move in (
                ["metropolized-gibbs"]
                if quick
                else ["gibbs", "metropolized-gibbs", "metropolis", "barker"]
            ):
                cases.append(
                    {
                        "name": f"{kind!s}-r{ranks:d}-{move!s}",
                        "kind": kind,
                        "ranks": ranks,
                        "move": move,
                    }
                )
    if not quick:
        if not skip_population:
            cases += [
                {
                    "name": "population",
                    "kind": "lambda",
                    "ranks": 1,
                    "move": "gibbs",
                    "population": True,
                }
            ]
        cases += [
            {
                "name": f"restart-{kind!s}-r{ranks:d}-at{at:d}",
                "kind": kind,
                "ranks": ranks,
                "move": "metropolized-gibbs",
                "restart": True,
                "restartAt": at,
            }
            for kind in ("metad", "metad-energy")
            for ranks in (1, 2)
            for at in (5, 7, 10)
        ]
        cases += [
            {
                "name": "energy-sparse",
                "kind": "energy",
                "ranks": 2,
                "move": "gibbs",
                "sparse": True,
            }
        ]
        cases += [
            {
                "name": kind + "-stride10",
                "kind": kind,
                "ranks": 2,
                "move": "gibbs",
                "stride": 10,
            }
            for kind in ("lambda", "energy")
        ]
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--provider")
    parser.add_argument("--backend", choices=["mpi", "tmpi"], default="mpi")
    parser.add_argument("--output", required=True)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--legacy", action="store_true", help="Use the legacy simulator"
    )
    parser.add_argument(
        "--skip-population",
        action="store_true",
        help="Skip the long CPU-only equilibrium control",
    )
    parser.add_argument("--filter", default="")
    parser.add_argument(
        "--launcher",
        default="/opt/local/bin/mpiexec-openmpi-mp --host localhost:{n} --bind-to none -n {n}"
        if sys.platform == "darwin"
        else "mpiexec -n {n}",
    )
    parser.add_argument("--mdrun-args", default="-nb cpu -pme cpu -update cpu")
    parser.add_argument("--serial-launcher", default="")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=180)
    options = parser.parse_args()
    os.environ["PLUMED_KERNEL"] = str(Path(options.kernel).resolve())
    os.environ["PLUMED_EXCHANGE_DEBUG"] = "1"
    if options.legacy:
        os.environ["GMX_DISABLE_MODULAR_SIMULATOR"] = "1"
    else:
        os.environ.pop("GMX_DISABLE_MODULAR_SIMULATOR", None)
    if options.backend == "tmpi":
        os.environ.pop("PLUMED_MPI_WRAPPER", None)
    if options.provider:
        os.environ["PLUMED_MPI_WRAPPER"] = str(Path(options.provider).resolve())
    args = SimpleNamespace(
        backend=options.backend,
        gmx=str(Path(options.gmx).resolve()),
        output=options.output,
        threads=options.threads,
        timeout=options.timeout,
        serial_launcher=options.serial_launcher,
        launcher=options.launcher,
        mdrun_args=options.mdrun_args,
        atol=3e-5,
        rtol=3e-5,
    )
    suite = module.Suite(args)
    results = []
    cases = cases_for(options.quick, options.skip_population)
    if options.legacy:
        # Native GROMACS disables legacy expanded-ensemble checkpoints (#4629).
        # Restart equivalence belongs to the supported modular-simulator suite.
        cases = [case for case in cases if not case.get("restart")]
    for case in cases:
        if options.filter and options.filter not in case["name"]:
            continue
        item = {"case": case}
        try:
            model = {
                "name": case["name"],
                "n": 1,
                "ranks": case["ranks"],
                "bias": "none",
                "lambdas": [0, 0.3, 0.7, 1],
            }
            root = suite.setup(model)
            directory = root / "rep0"
            nsteps = 20000 if case.get("population") else 31
            mdp = (
                (directory / "md.mdp")
                .read_text()
                .replace("integrator = md", "integrator = md-vv")
                .replace("free-energy = yes", "free-energy = expanded")
                .replace("nsteps = 31", f"nsteps = {nsteps:d}")
            )
            mdp += "\nnstexpanded = {:d}\nlmc-stats = no\nlmc-move = {!s}\nlmc-repeats = 3\nlmc-seed = 9831\nnsttcouple = 1\ninit-lambda-weights = 0 0.4 -0.2 0.1\n".format(
                1 if case.get("population") else 5, case["move"]
            )
            mdp = re.sub(r"(?m)^nstfout\s*=.*$", "nstfout = 1", mdp)
            if not re.search(r"(?m)^nstfout\s*=", mdp):
                mdp += "\nnstfout = 1\n"
            if case.get("sparse"):
                mdp = (
                    mdp.replace("nstcalcenergy = 1", "nstcalcenergy = 5")
                    .replace("nstenergy = 1", "nstenergy = 5")
                    .replace("nstlog = 1", "nstlog = 5")
                )
            (directory / "md.mdp").write_text(mdp)
            if case.get("population"):
                top = (
                    (directory / "topol.top")
                    .read_text()
                    .replace("0.34 2.0", "0.34 0.0")
                    .replace("0.34 8.0", "0.34 0.0")
                )
                (directory / "topol.top").write_text(top)
            kind = case["kind"]
            plumed = "CHECKPOINT FILE=plumed.cpt\nl: EXTRACV NAME=lambda\nd: DISTANCE ATOMS=1,2\n"
            fields = ["l", "d"]
            if "energy" in kind:
                plumed += "e: ENERGY\n"
                fields += ["e"]
            if kind in ("lambda", "lambda-energy"):
                plumed += "r: RESTRAINT ARG=l AT=1.3 KAPPA=2\n"
                fields += ["r.bias"]
            if kind in ("energy", "lambda-energy", "metad-energy"):
                plumed += "er: RESTRAINT ARG=e AT=0 KAPPA=2\n"
                fields += ["er.bias"]
            if kind.startswith("metad"):
                arguments = (
                    "l,d"
                    if kind == "metad2d"
                    else "l,e"
                    if kind == "metad-energy"
                    else "l"
                )
                sigma = (
                    "0.6,0.1"
                    if kind == "metad2d"
                    else "0.6,0.2"
                    if kind == "metad-energy"
                    else "0.6"
                )
                plumed += f"m: METAD ARG={arguments} SIGMA={sigma} HEIGHT=0.2 PACE=5 FILE=HILLS\n"
                fields += ["m.bias"]
            if case.get("stride"):
                plumed = re.sub(
                    r"((?:er|r): RESTRAINT[^\n]+)",
                    "\\1 STRIDE={:d}".format(case["stride"]),
                    plumed,
                )
            plumed += "PRINT ARG={} FILE=COLVAR STRIDE=1 FMT=%.12f\n".format(
                ",".join(fields)
            )
            (directory / "plumed.dat").write_text(plumed)
            suite.run(
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
                directory,
                "grompp-expanded.log",
            )
            command = (
                suite.launch(case["ranks"])
                + [
                    "-s",
                    "topol.tpr",
                    "-deffnm",
                    "run",
                    "-plumed",
                    "plumed.dat",
                    "-pin",
                    "off",
                    "-notunepme",
                ]
                + shlex.split(args.mdrun_args)
            )
            if case.get("restart"):
                suite.run(
                    command + ["-nsteps", str(case.get("restartAt", 10))],
                    directory,
                    "first.log",
                )
                suite.run(
                    command + ["-cpi", "run.cpt", "-append"], directory, "continued.log"
                )
            else:
                suite.run(command, directory, "mdrun.log")
            rows = module.analytic.table(next(directory.glob("COLVAR*")))
            if case.get("restart"):
                continuous = root / "continuous"
                continuous.mkdir()
                for name in ("topol.tpr", "plumed.dat"):
                    shutil.copyfile(directory / name, continuous / name)
                suite.run(command, continuous, "mdrun.log")
                reference = module.analytic.table(next(continuous.glob("COLVAR*")))
                if len(rows) != len(reference):
                    raise AssertionError("Restart/continuous output counts differ")
                maximum = 0.0
                for actual, expected in zip(rows, reference):
                    for field in actual:
                        error = abs(actual[field] - expected[field])
                        maximum = max(maximum, error)
                        if error > 3e-5 + 3e-5 * abs(expected[field]):
                            raise AssertionError(
                                "Restart differs from continuous run in " + field
                            )
                item["restart_max_error"] = maximum
            bystep = {round(row["time"] / 0.0001): row for row in rows}
            if sorted(bystep) != list(range(nsteps + 1)) or len(rows) != nsteps + 1:
                raise AssertionError("Missing/duplicate ordinary PLUMED output")
            for row in rows:
                if "e" in row:
                    expected = suite.physical(model, round(row["l"]), row["d"])
                    if abs(row["e"] - expected) > 3e-5 + 3e-5 * abs(expected):
                        raise AssertionError(
                            "Normal ENERGY is missing/stale on an MD step"
                        )
            if kind in ("energy", "lambda-energy"):
                suite.run([args.gmx, "dump", "-f", "run.trr"], directory, "forces.dump")
                forces = [
                    float(match[1].split(",")[0])
                    for match in re.finditer(
                        r"\bf\[\s*0\]\s*=\{([^}]+)\}",
                        (directory / "forces.dump").read_text(),
                    )
                ]
                if len(forces) != len(rows):
                    raise AssertionError("Missing force frames")
                maximum_force_error = 0.0
                for actual, row in zip(forces, rows):
                    state = round(row["l"])
                    distance = row["d"]
                    z = (0.34 / distance) ** 6
                    epsilon = 2 + 6 * model["lambdas"][state]
                    derivative = 24 * epsilon / distance * (z - 2 * z**2)
                    step = round(row["time"] / 0.0001)
                    stride = case.get("stride", 1)
                    expected = derivative * (
                        1 + (2 * row["e"] * stride if step % stride == 0 else 0)
                    )
                    error = abs(actual - expected)
                    maximum_force_error = max(maximum_force_error, error)
                    if error > 0.002 + 3e-4 * abs(expected):
                        raise AssertionError(
                            "Normal ENERGY force is wrong or an MCMC trial changed forces"
                        )
                item["force_max_error"] = maximum_force_error
            hills = (
                module.analytic.table(next(directory.glob("HILLS*")))
                if kind.startswith("metad")
                else []
            )
            if kind.startswith("metad") and len(hills) != nsteps // 5:
                raise AssertionError("Trial evaluation updated metadynamics history")

            def bias(state, energy, distance, step, kind=kind, hills=hills):
                result = (
                    (state - 1.3) ** 2 if kind in ("lambda", "lambda-energy") else 0
                )
                if kind in ("energy", "lambda-energy", "metad-energy"):
                    result += energy**2
                for hill in hills:
                    if round(hill["time"] / 0.0001) <= step:
                        squared = ((state - hill["l"]) / hill["sigma_l"]) ** 2
                        if kind == "metad2d":
                            squared += ((distance - hill["d"]) / hill["sigma_d"]) ** 2
                        if kind == "metad-energy":
                            squared += ((energy - hill["e"]) / hill["sigma_e"]) ** 2
                        result += hill["height"] * max(
                            0,
                            (math.exp(-0.5 * squared) - math.exp(-6.25))
                            / (1 - math.exp(-6.25)),
                        )
                return result

            checks = []
            log = (directory / "run.log").read_text()
            for match in re.finditer(
                r"PLUMED_LAMBDA step=(\d+) current=(\d+) state=(\d+) potential=(\S+) bias=(\S+) reduced=(\S+)",
                log,
            ):
                step, current, state = map(int, match.group(1, 2, 3))
                actual_u, actual_b, actual_reduced = map(float, match.group(4, 5, 6))
                row = bystep[step]
                if round(row["l"]) != current:
                    raise AssertionError("Normal lambda does not match the MD state")

                def potential(index, case=case, model=model, row=row):
                    return (
                        0
                        if case.get("population")
                        else suite.physical(model, index, row["d"])
                    )

                u, own = potential(state), potential(current)
                b = bias(state, u, row["d"], step)
                old_b = bias(current, own, row["d"], step)
                expected = (u - own + b - old_b) / (module.analytic.R * 300)
                error = max(
                    abs(actual_u - u), abs(actual_b - b), abs(actual_reduced - expected)
                )
                tolerance = 3e-5 + 3e-5 * max(abs(u), abs(b), abs(expected))
                if error > tolerance:
                    raise AssertionError(
                        f"Trial mismatch step{step:d} state{state:d}: error={error:g} tolerance={tolerance:g} U={actual_u:g}/{u:g} bias={actual_b:g}/{b:g} reduced={actual_reduced:g}/{expected:g}"
                    )
                checks.append(
                    {
                        "step": step,
                        "current": current,
                        "state": state,
                        "error": error,
                        "tolerance": tolerance,
                    }
                )
            if len(checks) < 4 * (nsteps // (1 if case.get("population") else 5)):
                raise AssertionError("Missing MCMC trial evaluations")
            for match in re.finditer(
                r"PLUMED_LAMBDA_MOVE step=(\d+) old=(\d+) new=(\d+)", log
            ):
                step, _old, new = map(int, match.group(1, 2, 3))
                if step + 1 in bystep and round(bystep[step + 1]["l"]) != new:
                    raise AssertionError(
                        "MCMC result was not applied on the next MD step"
                    )
            if case.get("population"):
                samples = rows[100:]
                counts = [
                    sum(round(row["l"]) == index for row in samples)
                    for index in range(4)
                ]
                weights = [
                    math.exp(offset - (index - 1.3) ** 2 / (module.analytic.R * 300))
                    for index, offset in enumerate([0, 0.4, -0.2, 0.1])
                ]
                expected = [value / sum(weights) for value in weights]
                for index, probability in enumerate(expected):
                    frequency = counts[index] / len(samples)
                    tolerance = (
                        6 * math.sqrt(probability * (1 - probability) / len(samples))
                        + 0.005
                    )
                    if abs(frequency - probability) > tolerance:
                        raise AssertionError("Wrong equilibrium lambda population")
                item["populations"] = {"counts": counts, "expected": expected}
            item.update(
                status="pass",
                trials=len(checks),
                max_error=max(check["error"] for check in checks),
            )
        except (
            AssertionError,
            RuntimeError,
            OSError,
            ValueError,
            KeyError,
            IndexError,
            subprocess.SubprocessError,
        ) as error:
            item.update(status="fail", error=str(error))
        results.append(item)
        (suite.root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
        print(case["name"], item["status"], item.get("error", ""), flush=True)
    return int(not results or any(item["status"] != "pass" for item in results))


if __name__ == "__main__":
    sys.exit(main())
