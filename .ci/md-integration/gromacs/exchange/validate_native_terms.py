"""Native alchemical acceptance retains mass/constraint deltaH terms with PLUMED."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import sys
from types import SimpleNamespace

location = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "native_lambda", location / "validate_lambda.py"
)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--launcher", required=True)
    options = parser.parse_args()
    os.environ["PLUMED_KERNEL"] = str(Path(options.kernel).resolve())
    args = SimpleNamespace(
        gmx=str(Path(options.gmx).resolve()),
        output=options.output,
        threads=1,
        timeout=120,
        serial_launcher="",
        launcher=options.launcher,
        mdrun_args="-nb cpu -pme cpu -update cpu",
        atol=3e-5,
        rtol=3e-5,
    )
    suite = native.Suite(args)
    results = []
    for constraints in [False, True]:
        case = dict(
            name=f"native-mass-constraints{int(constraints)}",
            n=2,
            ranks=1,
            bias="none",
            nstdhdl=1,
        )
        root = suite.setup(case)
        for replica in range(2):
            directory = root / f"rep{replica}"
            top = (
                (directory / "topol.top")
                .read_text()
                .replace("ArB 0.0 39.948", "ArB 0.0 59.922")
            )
            if constraints:
                top = top.replace("AR 1\n[ atoms ]", "AR 0\n[ atoms ]")
                top = top.replace(
                    "[ system ]",
                    "2 ArA 1 AR AR 2 0.0 39.948 ArB 0.0 59.922\n[ constraints ]\n1 2 1 0.45 0.50\n[ system ]",
                )
                top = top.replace("AR 2\n", "AR 1\n")
            (directory / "topol.top").write_text(top)
            mdp = (
                directory / "md.mdp"
            ).read_text() + "\ndhdl-derivatives = no\ndhdl-print-energy = potential\n"
            (directory / "md.mdp").write_text(mdp)
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
                "grompp-native.log",
            )
        common = (
            suite.launch(2)
            + [
                "-multidir",
                "rep0",
                "rep1",
                "-s",
                "topol.tpr",
                "-replex",
                "5",
                "-reseed",
                "2026",
                "-notunepme",
                "-pin",
                "off",
            ]
            + shlex.split(args.mdrun_args)
        )
        suite.run(common + ["-deffnm", "plain"], root, "plain.log")
        suite.run(
            common + ["-deffnm", "biased", "-plumed", "plumed.dat"], root, "biased.log"
        )
        # Native dhdl columns give the full deltaH, including mass/constraint terms.
        tables = []
        for replica in range(2):
            file = root / f"rep{replica}/biased.xvg"
            data = [
                line.split()
                for line in file.read_text().splitlines()
                if line and not line.startswith(("#", "@"))
            ]
            tables.append(
                {round(float(row[0]) / 0.0001): list(map(float, row)) for row in data}
            )
            plain = (root / f"rep{replica}/plain.log").read_text()
            biased = (root / f"rep{replica}/biased.log").read_text()
            for prefix in ["Repl ex", "Repl pr"]:
                if re.findall(rf"^{prefix}.*$", plain, re.M) != re.findall(
                    rf"^{prefix}.*$", biased, re.M
                ):
                    raise AssertionError(
                        "Zero PLUMED bias changed native exchange decisions/probabilities"
                    )
        checks = []
        for step, correction, delta in re.findall(
            r"PLUMED_GREX step=(\d+) a=0 b=1 correction=(\S+) delta=(\S+)",
            (root / "rep0/biased.log").read_text(),
        ):
            step = int(step)
            # Last two entries are the two foreign-state deltaH values.
            left, right = tables[0][step], tables[1][step]
            expected = (left[-1] - left[-2] + right[-2] - right[-1]) / (
                native.analytic.R * 300
            )
            error = abs(float(delta) - expected)
            if abs(float(correction)) > 1e-6 or error > 2e-4 + 1e-4 * abs(expected):
                raise AssertionError(
                    f"Native deltaH terms lost: step{step} delta{delta} expected{expected}"
                )
            checks.append(
                dict(step=step, delta=float(delta), expected=expected, error=error)
            )
        if not checks:
            raise AssertionError("No native acceptance checks")
        results.append(dict(case=case, status="pass", checks=checks))
        (Path(args.output) / "results.json").write_text(
            json.dumps(results, indent=2) + "\n"
        )
        print(case["name"], "pass", flush=True)
    # This combination is not supported by the native acceptance machinery.
    # The PLUMED integration must reject it instead of appearing to implement it.
    case = dict(name="expanded-replex-rejected", n=2, ranks=1, bias="none")
    root = suite.setup(case)
    for replica in range(2):
        directory = root / f"rep{replica}"
        path = directory / "md.mdp"
        text = path.read_text().replace("integrator = md\n", "integrator = md-vv\n")
        text = text.replace("free-energy = yes", "free-energy = expanded")
        path.write_text(text + "\nnstexpanded = 1\ninit-lambda-weights = 0 0\n")
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
    code = suite.run(
        common + ["-deffnm", "run", "-plumed", "plumed.dat"],
        root,
        "rejected.log",
        check=False,
    )
    diagnostic = " ".join((root / "rejected.log").read_text().split())
    if (
        not code
        or "Expanded-ensemble lambda MCMC combined with replica exchange is unsupported"
        not in diagnostic
    ):
        raise AssertionError(
            "Unsupported MCMC plus replica exchange did not fail clearly"
        )
    results.append(
        dict(case=case, status="pass", expected_failure=True, exit_code=code)
    )
    (Path(args.output) / "results.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(case["name"], "pass", flush=True)
    # Preserve native startup validation when the PLUMED provider is inactive.
    case = dict(name="identical-native-replicas-rejected", n=2, ranks=1, bias="none")
    root = native.analytic.Suite.setup(suite, case)
    code = suite.run(
        common + ["-deffnm", "run"], root, "rejected.log", check=False
    )
    diagnostic = " ".join((root / "rejected.log").read_text().split()).lower()
    if not code or not any(
        text in diagnostic for text in ("nothing to exchange", "identical")
    ):
        raise AssertionError("Identical native replicas were accepted without PLUMED")
    results.append(
        dict(case=case, status="pass", expected_failure=True, exit_code=code)
    )
    (Path(args.output) / "results.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(case["name"], "pass", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
