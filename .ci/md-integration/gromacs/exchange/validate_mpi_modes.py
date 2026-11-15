"""Check thread-MPI PLUMED MCMC and ordinary-run controls."""

import argparse
import importlib.util
import json
import math
import os
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location(
    "lambda_suite", Path(__file__).with_name("validate_lambda.py")
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True, help="GROMACS built with thread-MPI")
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    os.environ["PLUMED_KERNEL"] = str(Path(args.kernel).resolve())
    os.environ.pop("PLUMED_MPI_WRAPPER", None)
    options = SimpleNamespace(
        gmx=str(Path(args.gmx).resolve()),
        output=args.output,
        threads=1,
        timeout=120,
        serial_launcher="",
        launcher="",
    )
    suite = module.Suite(options)
    suite.env.pop("GMX_DISABLE_MODULAR_SIMULATOR", None)
    template = (
        suite.setup(
            {
                "name": "inputs",
                "n": 1,
                "ranks": 1,
                "bias": "none",
                "lambdas": [0, 0.3, 0.7, 1],
            }
        )
        / "rep0"
    )
    base_mdp = (template / "md.mdp").read_text()
    results = []

    def run_case(
        name, ranks, integrator="md", bias="coordinate", expanded=False, legacy=False
    ):
        item = {
            "name": name,
            "ranks": ranks,
            "integrator": integrator,
            "bias": bias,
            "expanded": expanded,
            "legacy": legacy,
        }
        try:
            directory = suite.root / name
            directory.mkdir()
            for filename in ("conf.gro", "topol.top"):
                shutil.copyfile(template / filename, directory / filename)
            mdp = base_mdp.replace("integrator = md", "integrator = " + integrator)
            if expanded:
                mdp = mdp.replace("free-energy = yes", "free-energy = expanded")
                mdp += "\nnstexpanded = 5\nlmc-stats = no\nlmc-move = metropolized-gibbs\nlmc-seed = 9831\n"
            (directory / "md.mdp").write_text(mdp)
            coordinate = "d: DISTANCE ATOMS=1,2\nr: RESTRAINT ARG=d AT=0.4 KAPPA=2\nPRINT ARG=d,r.bias FILE=COLVAR STRIDE=1 FMT=%.12f\n"
            energy = "d: DISTANCE ATOMS=1,2\ne: ENERGY\nr: RESTRAINT ARG=e AT=0 KAPPA=2\nPRINT ARG=d,e,r.bias FILE=COLVAR STRIDE=1 FMT=%.12f\n"
            (directory / "plumed.dat").write_text(
                energy if bias == "energy" else coordinate
            )
            suite.run(
                [
                    options.gmx,
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
                "grompp.log",
            )
            command = [
                options.gmx,
                "mdrun",
                "-s",
                "topol.tpr",
                "-deffnm",
                "run",
                "-ntmpi",
                str(ranks),
                "-ntomp",
                "1",
                "-pin",
                "off",
                "-notunepme",
                "-nb",
                "cpu",
                "-pme",
                "cpu",
                "-update",
                "cpu",
            ]
            if bias != "none":
                command += ["-plumed", "plumed.dat"]
            if legacy:
                suite.env["GMX_DISABLE_MODULAR_SIMULATOR"] = "1"
            else:
                suite.env.pop("GMX_DISABLE_MODULAR_SIMULATOR", None)
            rc = suite.run(command, directory, "mdrun.log", check=False)
            item["exit_code"] = rc
            log = (directory / "mdrun.log").read_text()
            assert rc == 0, log[-4000:]
            assert (directory / "run.gro").is_file(), "Missing final configuration"
            if bias != "none":
                rows = module.analytic.table(directory / "COLVAR")
                assert len(rows) == 32, "Missing ordinary PLUMED steps"
                for row in rows:
                    expected = (
                        (row["d"] - 0.4) ** 2 if bias == "coordinate" else row["e"] ** 2
                    )
                    assert math.isclose(
                        row["r.bias"], expected, abs_tol=3e-5, rel_tol=3e-5
                    ), row
                item["rows"] = len(rows)
            item["status"] = "pass"
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
        print(name, item["status"], item.get("error", ""), flush=True)

    for ranks in (1, 2):
        for legacy in (False, True):
            run_case(
                "biased-mcmc-r{:d}-{!s}".format(
                    ranks, "legacy" if legacy else "modular"
                ),
                ranks,
                integrator="md-vv",
                bias="energy",
                expanded=True,
                legacy=legacy,
            )
        run_case(
            f"native-mcmc-r{ranks:d}",
            ranks,
            integrator="md-vv",
            bias="none",
            expanded=True,
        )
        for integrator in ("md", "md-vv"):
            for bias in ("coordinate", "energy"):
                run_case(
                    f"ordinary-{integrator!s}-{bias!s}-r{ranks:d}",
                    ranks,
                    integrator,
                    bias,
                )
    return int(any(item["status"] != "pass" for item in results))


if __name__ == "__main__":
    raise SystemExit(main())
