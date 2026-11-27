"""Check selective local-auxiliary guards and supported independent/walker MD."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--launcher", required=True)
    options = parser.parse_args()
    gmx = str(Path(options.gmx).resolve())
    root = Path(options.output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    location = Path(__file__).resolve().parent
    model = location / "multisim-auxiliary"
    env = os.environ | {
        "PLUMED_KERNEL": str(Path(options.kernel).resolve()),
        "OMP_NUM_THREADS": "1",
        "PLUMED_NUM_THREADS": "1",
    }

    def run(command, cwd, filename):
        with (cwd / filename).open("w") as log:
            process = subprocess.Popen(
                command,
                cwd=cwd,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            try:
                return process.wait(timeout=120)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait(timeout=5)
                raise AssertionError(
                    f"{filename}: coupled evaluation did not complete"
                ) from None

    results = []
    for case in [
        "ordinary",
        "mcmc",
        "shells",
        "mcmc-independent",
        "shells-independent",
        "mcmc-walkers",
        "shells-walkers",
    ]:
        work = root / case
        work.mkdir(exist_ok=True)
        for replica in range(2):
            rep = work / f"rep{replica}"
            rep.mkdir(exist_ok=True)
            template = location / "shell" if case.startswith("shells") else model
            shutil.copy2(template / "topol.top", rep)
            shutil.copy2(template / "conf.gro", rep)
            mdp = (template / "md.mdp").read_text()
            if case == "ordinary":
                ignored = [
                    "free-energy",
                    "init-lambda",
                    "fep-lambdas",
                    "calc-lambda",
                    "nstdhdl",
                    "sc-",
                    "nstexpanded",
                    "lmc-",
                    "init-lambda-weights",
                ]
                mdp = (
                    "\n".join(
                        line
                        for line in mdp.splitlines()
                        if not any(line.startswith(key) for key in ignored)
                    )
                    + "\n"
                )
                mdp = mdp.replace("nsteps = 6", "nsteps = 1")
                plumed = "d: DISTANCE ATOMS=1,2 NOPBC\ne: ENSEMBLE ARG=d\nr: RESTRAINT ARG=e.d AT=0 KAPPA=2\nPRINT ARG=d,e.d,r.bias FILE=COLVAR\n"
            elif case.startswith("mcmc"):
                mdp = mdp.replace(
                    "init-lambda-state = 0", f"init-lambda-state = {replica}"
                )
                plumed = "l: EXTRACV NAME=lambda\ne: ENSEMBLE ARG=l\nr: RESTRAINT ARG=e.l AT=0 KAPPA=2\n"
            else:
                # One replica converges at its first force call; the other needs
                # more. An unguarded ensemble query then mismatches collectives.
                if case == "shells":
                    mdp = mdp.replace("nsteps = 2", "nsteps = 0")
                xyz = "Shell collective control\n2\n"
                for atom, atomname, x in [
                    (1, "C", 0),
                    (2, "S", 1e-7 if replica == 0 else 0.1),
                ]:
                    xyz += f"{1:5d}{'POL':<5}{atomname:>5}{atom:5d}{x:13.8f}{1.0:13.8f}{1.0:13.8f}\n"
                (rep / "conf.gro").write_text(xyz + "3.0 3.0 3.0\n")
                plumed = "d: DISTANCE ATOMS=1,2 NOPBC\ne: ENSEMBLE ARG=d\nr: RESTRAINT ARG=e.d AT=0 KAPPA=0\n"
            if case.endswith("-independent") or case.endswith("-walkers"):
                if case.startswith("mcmc"):
                    variable, definition = "l", "l: EXTRACV NAME=lambda\n"
                else:
                    variable, definition = "d", "d: DISTANCE ATOMS=1,2 NOPBC\n"
                if case.endswith("-walkers"):
                    bias = f"m: METAD ARG={variable} SIGMA=1 HEIGHT=0.01 PACE=1 WALKERS_MPI FILE=HILLS\n"
                else:
                    bias = f"m: RESTRAINT ARG={variable} AT=0 KAPPA=2\n"
                plumed = (
                    definition + bias + f"PRINT ARG={variable},m.bias FILE=COLVAR\n"
                )
            (rep / "md.mdp").write_text(mdp)
            (rep / "plumed.dat").write_text(plumed)
            if run(
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
                rep,
                "grompp.log",
            ):
                raise AssertionError(f"{case}: invalid input; inspect grompp.log")
        command = shlex.split(options.launcher.format(n=2)) + [
            gmx,
            "mdrun",
            "-multidir",
            "rep0",
            "rep1",
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
        code = run(command, work, "mdrun.log")
        log = (work / "mdrun.log").read_text()
        if case == "ordinary":
            if code != 0:
                raise AssertionError(
                    "ordinary multi-simulation ensemble MD was rejected"
                )
            for replica in range(2):
                rep = work / f"rep{replica}"
                if "averaging over 2 replicas" not in (rep / "PLUMED.OUT").read_text():
                    raise AssertionError(
                        "ordinary control did not exercise replica coupling"
                    )
                files = list(rep.glob("COLVAR*"))
                rows = [
                    line
                    for line in files[0].read_text().splitlines()
                    if line and not line.startswith("#")
                ]
                if len(rows) != 2:
                    raise AssertionError("ordinary control did not finish both steps")
        elif case in ["mcmc", "shells"]:
            required = "action e (ENSEMBLE) communicates between replicas"
            mode = 2 if case == "mcmc" else 1
            normalized = " ".join(log.split())
            if (
                code == 0
                or required not in normalized
                or f"local auxiliary mode {mode}" not in normalized
            ):
                raise AssertionError(f"{case}: missing action-specific rejection")
        else:
            if code != 0:
                raise AssertionError(
                    f"{case}: supported local auxiliary evaluation failed"
                )
            for replica in range(2):
                rep = work / f"rep{replica}"
                files = list(rep.glob("COLVAR*"))
                rows = (
                    [
                        line
                        for line in files[0].read_text().splitlines()
                        if line and not line.startswith("#")
                    ]
                    if files
                    else []
                )
                expected = 7 if case.startswith("mcmc") else 3
                if len(rows) != expected:
                    raise AssertionError(f"{case}: missing ordinary output records")
            if case.endswith("-walkers"):
                # WALKERS_MPI writes one shared HILLS file from the leader.
                # Count both walkers' deposits, excluding the initial step and
                # every local auxiliary force/lambda evaluation.
                hills = list(work.glob("rep*/HILLS*"))
                deposits = sum(
                    sum(
                        line and not line.startswith("#")
                        for line in path.read_text().splitlines()
                    )
                    for path in hills
                )
                if deposits != 2 * (expected - 1):
                    raise AssertionError(
                        f"{case}: wrong number of committed walker hills"
                    )
        results.append({"case": case, "status": "pass", "returncode": code})
        print(f"{case}: PASS", flush=True)
    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
