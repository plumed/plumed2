"""Native lambda exchange: independent LJ/softcore and PLUMED bias oracles."""

import argparse
import importlib.util
import json
import math
import os
import re
import sys
from pathlib import Path
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location(
    "analytic", Path(__file__).with_name("validate.py")
)
analytic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analytic)


class Suite(analytic.Suite):
    def setup(self, case):
        root = super().setup(case)
        lambdas = (
            case["lambdas"]
            if "lambdas" in case
            else [i / (case["n"] - 1) for i in range(case["n"])]
        )
        for i in range(case["n"]):
            directory = root / (f"rep{i:d}")
            top = "[ defaults ]\n1 2 no\n[ atomtypes ]\nArA 39.948 0.0 A 0.34 2.0\nArB 39.948 0.0 A 0.34 {!s}\n[ moleculetype ]\nAR 1\n[ atoms ]\n1 ArA 1 AR AR 1 0.0 39.948 ArB 0.0 39.948\n[ system ]\nAlchemical LJ pair\n[ molecules ]\nAR 2\n".format(
                0.0 if case.get("softcore") else 8.0
            )
            (directory / "topol.top").write_text(top)
            schedule = " ".join(str(value) for value in lambdas)
            mdp = (directory / "md.mdp").read_text()
            mdp += "\nfree-energy = yes\ninit-lambda-state = {:d}\nfep-lambdas = {!s}\ncalc-lambda-neighbors = -1\nnstdhdl = {:d}\nsc-alpha = {!s}\nsc-power = 1\nsc-sigma = 0.3\n".format(
                i, schedule, case.get("nstdhdl", 10), 0.5 if case.get("softcore") else 0
            )
            if "lambdaNeighbors" in case:
                mdp = mdp.replace(
                    "calc-lambda-neighbors = -1",
                    "calc-lambda-neighbors = {:d}".format(case["lambdaNeighbors"]),
                )
            if case.get("dispCorr"):
                mdp += "\nDispCorr = EnerPres\n"
            (directory / "md.mdp").write_text(mdp)
            if case.get("pme"):
                top = top.replace("AR 1\n[ atoms ]", "AR 0\n[ atoms ]")
                top = top.replace(
                    "1 ArA 1 AR AR 1 0.0 39.948 ArB 0.0 39.948",
                    "1 ArA 1 AR AR 1 0.04 39.948 ArB 0.08 39.948\n2 ArA 1 AR AR 1 -0.04 39.948 ArB -0.08 39.948",
                )
                top = top.replace("AR 2\n", "AR 1\n")
                (directory / "topol.top").write_text(top)
                mdp = mdp.replace(
                    "coulombtype = Cut-off",
                    "coulombtype = PME\nfourier-nx = 32\nfourier-ny = 32\nfourier-nz = 32\npme-order = 4\newald-rtol = 1e-6",
                )
                (directory / "md.mdp").write_text(mdp)
            # Keep ENERGY forces moderate while making wrong foreign values detectable.
            plumed = (
                (directory / "plumed.dat")
                .read_text()
                .replace("KAPPA=100\n", "KAPPA=2\n")
            )
            plumed = plumed.replace("SIGMA=0.1", "SIGMA=0.2")
            if case["bias"] in ("lambda-coordinate", "lambda-energy"):
                plumed = f"l: CONSTANT VALUE={i}\nd: DISTANCE ATOMS=1,2\n"
                if case["bias"] == "lambda-coordinate":
                    plumed += "f: CUSTOM ARG=l,d FUNC=x*y-0.1 PERIODIC=NO\n"
                else:
                    plumed += "e: ENERGY\nf: CUSTOM ARG=e,l FUNC=x-0.3*y PERIODIC=NO\n"
                plumed += "r: RESTRAINT ARG=f AT=0 KAPPA=2\nPRINT ARG=d,l,r.bias FILE=COLVAR STRIDE=1 FMT=%.12f\n"
            if case.get("stride"):
                plumed = re.sub(
                    r"((?:er|r): RESTRAINT[^\n]+)",
                    "\\1 STRIDE={:d}".format(case["stride"]),
                    plumed,
                )
                plumed = re.sub(r"PRINT ARG=\S+", "PRINT ARG=d", plumed)
            (directory / "plumed.dat").write_text(plumed)
            self.run(
                [
                    self.args.gmx,
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
                "grompp-lambda.log",
            )
        return root

    def physical(self, case, replica, distance):
        lam = (
            case["lambdas"]
            if "lambdas" in case
            else [i / (case["n"] - 1) for i in range(case["n"])]
        )[replica]
        sigma6 = 0.34**6
        if case.get("softcore"):
            z = sigma6 / (distance**6 + 0.5 * lam * sigma6)
            return 8 * (1 - lam) * (z * z - z)
        z = sigma6 / distance**6
        return 4 * (2 + 6 * lam) * (z * z - z)

    def bias(self, case, replica, distance, step, hills, energy=None):
        kind = case["bias"]
        energy = self.physical(case, replica, distance) if energy is None else energy
        result = 0.0
        if kind == "lambda-coordinate":
            result += (replica * distance - 0.1) ** 2
        if kind == "lambda-energy":
            result += (energy - 0.3 * replica) ** 2
        if kind in ("coordinate", "combined"):
            result += (distance - (0.3 + 0.2 * replica)) ** 2
        if kind in ("energy", "combined", "metad-energy"):
            result += (energy - 0.05 * replica) ** 2
        if kind in ("metad", "metad-energy"):
            argument = "e" if kind == "metad-energy" else "d"
            value = energy if argument == "e" else distance
            for hill in hills[replica]:
                hill_step = round(hill["time"] / 0.0001)
                if hill_step < step or (
                    getattr(self.args, "bias_order", "pre") == "post"
                    and hill_step == step
                ):
                    scaled = (value - hill[argument]) / hill["sigma_" + argument]
                    result += hill["height"] * max(
                        0.0,
                        (math.exp(-0.5 * scaled**2) - math.exp(-6.25))
                        / (1 - math.exp(-6.25)),
                    )
        return result

    def verify(self, case, root, prefix="run"):
        rows = [
            analytic.table(next((root / (f"rep{i:d}")).glob("COLVAR*")))
            for i in range(case["n"])
        ]
        bystep = [
            {round(row["time"] / 0.0001): row for row in replica} for replica in rows
        ]
        hills = [
            analytic.table(next((root / (f"rep{i:d}")).glob("HILLS*")))
            if case["bias"] in ("metad", "metad-energy")
            else []
            for i in range(case["n"])
        ]
        cross = {}
        if case.get("pme"):
            for i in range(case["n"]):
                directory = root / (f"rep{i:d}")
                for j in range(case["n"]):
                    name = f"oracle-{j:d}"
                    command = self.launch(1) + [
                        "-s",
                        "topol.tpr",
                        "-rerun",
                        str(root / (f"rep{j:d}") / (prefix + ".trr")),
                        "-deffnm",
                        name,
                        "-nb",
                        "cpu",
                        "-pme",
                        "cpu",
                        "-update",
                        "cpu",
                        "-notunepme",
                        "-pin",
                        "off",
                    ]
                    self.run(command, directory, name + ".log")
                    self.run(
                        [
                            self.args.gmx,
                            "energy",
                            "-f",
                            name + ".edr",
                            "-o",
                            name + ".xvg",
                            "-xvg",
                            "none",
                            "-dp",
                        ],
                        directory,
                        name + "-energy.log",
                        stdin="Potential\n0\n",
                    )
                    cross[i, j] = {
                        round(float(line.split()[0]) / 0.0001): float(line.split()[1])
                        for line in (directory / (name + ".xvg"))
                        .read_text()
                        .splitlines()
                        if line and not line.startswith(("#", "@"))
                    }
        events = {}
        for log in root.glob("rep*/" + prefix + ".log"):
            step = pair = None
            for line in log.read_text().splitlines():
                match = re.search(r"Replica exchange at step (\d+)", line)
                if match:
                    step, pair = int(match[1]), None
                match = re.match(
                    r"Repl\s+(\d+)\s*<->\s*(\d+).*dE_term\s*=\s*(\S+)", line
                )
                if match:
                    pair = int(match[1]), int(match[2])
                match = re.search(r"dplumed\s*=\s*(\S+)\s+dE_Term\s*=\s*(\S+)", line)
                if match and step is not None and pair:
                    rounded = float(match[2])
                    digits = int(match[2].lower().split("e")[1])
                    events[(step,) + pair] = (rounded, 0.5 * 10 ** (digits - 3))
                match = re.search(
                    r"PLUMED_GREX step=(\d+) a=(\d+) b=(\d+) correction=(\S+) delta=(\S+)",
                    line,
                )
                if match:
                    events[tuple(map(int, match.group(1, 2, 3)))] = (
                        float(match[5]),
                        0.0,
                    )
        if not events:
            raise AssertionError("No exchange diagnostics")
        checks = []
        for (step, i, j), (actual, rounding) in sorted(events.items()):
            xi, xj = bystep[i][step]["d"], bystep[j][step]["d"]

            def beta(k):
                return 1 / (
                    analytic.R * (300 + 100 * k if case.get("temperatures") else 300)
                )

            def u(h, configuration, step=step):
                return (
                    cross[h, configuration][step]
                    if cross
                    else self.physical(case, h, bystep[configuration][step]["d"])
                )

            def b(h, configuration, step=step, u=u):
                return self.bias(
                    case,
                    h,
                    bystep[configuration][step]["d"],
                    step,
                    hills,
                    u(h, configuration),
                )

            physical = beta(i) * (u(i, j) - u(i, i)) + beta(j) * (u(j, i) - u(j, j))
            bias = beta(i) * (b(i, j) - b(i, i)) + beta(j) * (b(j, i) - b(j, j))
            expected = physical + bias
            tolerance = self.args.atol + self.args.rtol * abs(expected) + rounding
            error = abs(actual - expected)
            for replica, distance in ((i, xi), (j, xj)):
                if (
                    "e" in bystep[replica][step]
                    and abs(bystep[replica][step]["e"] - u(replica, replica))
                    > tolerance
                ):
                    raise AssertionError(
                        "Normal ENERGY does not match the analytic Hamiltonian"
                    )
            checks.append(
                {
                    "step": step,
                    "i": i,
                    "j": j,
                    "physical": physical,
                    "bias": bias,
                    "actual": actual,
                    "expected": expected,
                    "error": error,
                    "tolerance": tolerance,
                }
            )
        failures = [check for check in checks if check["error"] > check["tolerance"]]
        if failures:
            raise AssertionError(json.dumps(failures[:3]))
        return checks


def cases_for(quick=False, pme=False):
    cases = []
    for kind in (
        "none",
        "coordinate",
        "energy",
        "metad",
        "metad-energy",
        "lambda-coordinate",
        "lambda-energy",
    ):
        for softcore in [False] if quick else [False, True]:
            for ranks in [1] if quick else [1, 2]:
                for temperatures in [False] if quick else [False, True]:
                    name = f"{kind!s}-sc{softcore:d}-r{ranks:d}-t{temperatures:d}"
                    cases.append(
                        {
                            "name": name,
                            "n": 2,
                            "ranks": ranks,
                            "bias": kind,
                            "softcore": softcore,
                            "temperatures": temperatures,
                        }
                    )
    if not quick:
        for kind in (
            "coordinate",
            "energy",
            "metad-energy",
            "lambda-coordinate",
            "lambda-energy",
        ):
            cases += [
                {"name": kind + "-three-replica", "n": 3, "ranks": 1, "bias": kind},
                {
                    "name": kind + "-four-scrambled",
                    "n": 4,
                    "ranks": 1,
                    "bias": kind,
                    "lambdas": [0.8, 0, 1, 0.3],
                    "custom": True,
                },
                {
                    "name": kind + "-restart",
                    "n": 2,
                    "ranks": 1,
                    "bias": kind,
                    "restart": True,
                },
                {
                    "name": kind + "-staggered",
                    "n": 2,
                    "ranks": 1,
                    "bias": kind,
                    "pace": 7,
                },
            ]
        cases.append(
            {
                "name": "missing-foreign-lambda-rejected",
                "n": 2,
                "ranks": 1,
                "bias": "energy",
                "lambdaNeighbors": 0,
                "expected_failure": "selected partner's lambda energy",
            }
        )
        for kind in ("energy", "metad-energy"):
            cases.append(
                {
                    "name": kind + "-neighbor-table",
                    "n": 3,
                    "ranks": 2,
                    "bias": kind,
                    "lambdaNeighbors": 1,
                }
            )
        for kind in ("energy", "lambda-coordinate", "lambda-energy"):
            cases.append(
                {
                    "name": kind + "-stride10",
                    "n": 2,
                    "ranks": 1,
                    "bias": kind,
                    "stride": 10,
                    "expected_failure": "requires bias force STRIDE to divide the exchange interval",
                }
            )
    if not quick and not pme:
        cases += [
            dict(name=kind + "-stride5", n=2, ranks=1, bias=kind, stride=5)
            for kind in ("energy", "lambda-coordinate", "lambda-energy")
        ]
        cases.append(
            dict(
                name="hrex-free-energy-rejected",
                n=2,
                ranks=1,
                bias="energy",
                hrex=True,
                expected_failure="-hrex cannot be combined with native free-energy mode",
            )
        )
    if pme:
        cases = [
            {
                "name": f"pme-{kind!s}-r{ranks:d}",
                "n": 2,
                "ranks": ranks,
                "bias": kind,
                "temperatures": True,
                "pme": True,
            }
            for kind in ("coordinate", "energy", "metad-energy")
            for ranks in (1, 2)
        ]
        cases += [
            {
                "name": f"pme-{kind}-tail-r2",
                "n": 2,
                "ranks": 2,
                "bias": kind,
                "temperatures": True,
                "pme": True,
                "dispCorr": True,
            }
            for kind in ("energy", "metad-energy")
        ]
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--provider")
    parser.add_argument("--output", required=True)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--pme", action="store_true")
    parser.add_argument("--expect-energy-failure", action="store_true")
    parser.add_argument("--filter", default="")
    parser.add_argument("--bias-order", choices=["pre", "post"], default="pre")
    parser.add_argument(
        "--launcher",
        default="/opt/local/bin/mpiexec-openmpi-mp --host localhost:{n} --bind-to none -n {n}"
        if sys.platform == "darwin"
        else "mpiexec -n {n}",
    )
    parser.add_argument("--mdrun-args", default="-nb cpu -pme cpu -update cpu")
    parser.add_argument("--serial-launcher", default="")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=120)
    options = parser.parse_args()
    os.environ["PLUMED_KERNEL"] = str(Path(options.kernel).resolve())
    if options.provider:
        os.environ["PLUMED_MPI_WRAPPER"] = str(Path(options.provider).resolve())
    args = SimpleNamespace(
        gmx=str(Path(options.gmx).resolve()),
        output=options.output,
        threads=options.threads,
        timeout=options.timeout,
        serial_launcher=options.serial_launcher,
        launcher=options.launcher,
        mdrun_args=options.mdrun_args,
        atol=3e-5,
        rtol=3e-5,
        bias_order=options.bias_order,
    )
    suite = Suite(args)
    cases = cases_for(options.quick, options.pme)
    results = []
    for case in cases:
        if options.filter and options.filter not in case["name"]:
            continue
        item = suite.case(case)
        if options.expect_energy_failure and case["bias"] in ("energy", "metad-energy"):
            item["expected_legacy_error"] = item["status"] == "fail" and item.get(
                "error", ""
            ).startswith('[{"step"')
        results.append(item)
        (suite.root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
        print(case["name"], item["status"], item.get("error", ""), flush=True)
    if not results:
        raise RuntimeError("No cases selected")
    return int(
        any(
            item["status"] != "pass" and not item.get("expected_legacy_error")
            for item in results
        )
    )


if __name__ == "__main__":
    sys.exit(main())
