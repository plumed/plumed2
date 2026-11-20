"""Portable GREX integration tests. No third-party Python dependencies.
Independent LJ/RESTRAINT/METAD energy oracles; never use the implementation's
reported correction as its own reference. Every failing case is retained.
"""

import argparse
import contextlib
import hashlib
import json
import math
import os
import platform
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

R = 0.00831446261815324


def table(path):
    fields = None
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#! FIELDS"):
            fields = line.split()[2:]
        elif line and not line.startswith("#"):
            out.append(dict(zip(fields, map(float, line.split()))))
    return out


def failure_diagnostics(root):
    """Keep kernel errors even when MPI truncates the launcher's stack dump."""
    logs = [
        root / "mdrun.log",
        *root.glob("rep*/run.log"),
        *root.glob("rep*/PLUMED.OUT*"),
    ]
    return "\n".join(path.read_text(errors="replace") for path in logs)


class Suite:
    def __init__(self, args):
        self.args = args
        self.root = Path(args.output).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.env = os.environ.copy()
        self.env.update(OMP_NUM_THREADS=str(args.threads), GMX_MAXBACKUP="-1")
        self.results = []
        self.run([self.args.gmx, "--version"], self.root, "gromacs-version.txt")
        kernel = self.env.get("PLUMED_KERNEL", "")
        metadata = {
            "arguments": vars(args),
            "python": sys.version,
            "platform": platform.platform(),
            "kernel": kernel,
        }
        if kernel and Path(kernel).is_file():
            metadata["kernel_sha256"] = hashlib.sha256(
                Path(kernel).read_bytes()
            ).hexdigest()
        (self.root / "environment.json").write_text(
            json.dumps(metadata, indent=2) + "\n"
        )

    def run(self, cmd, cwd, log, stdin=None, check=True, separate_stderr=False):
        if cmd[0] == self.args.gmx and getattr(self.args, "serial_launcher", ""):
            cmd = shlex.split(self.args.serial_launcher) + cmd
        with (
            (cwd / log).open("w") as out,
            (
                (cwd / (log + ".stderr")).open("w")
                if separate_stderr
                else contextlib.nullcontext(subprocess.STDOUT)
            ) as err,
        ):
            r = subprocess.run(
                cmd,
                cwd=cwd,
                env=self.env,
                input=stdin,
                text=True,
                stdout=out,
                stderr=err,
                timeout=self.args.timeout,
                check=False,
            )
        if check and r.returncode:
            raise RuntimeError(f"{log}: exit {r.returncode}; inspect {cwd / log}")
        return r.returncode

    def launch(self, n):
        if getattr(self.args, "backend", "mpi") == "tmpi":
            return [
                self.args.gmx,
                "mdrun",
                "-ntmpi",
                str(n),
                "-ntomp",
                str(self.args.threads),
            ]
        return [x.format(n=n) for x in shlex.split(self.args.launcher)] + [
            self.args.gmx,
            "mdrun",
        ]

    def setup(self, c):
        root = self.root / c["name"]
        if root.exists():
            raise FileExistsError(
                f"Case directory already exists: {root}; use a fresh output directory"
            )
        root.mkdir()
        for i in range(c["n"]):
            d = root / f"rep{i}"
            d.mkdir(exist_ok=True)
            t = 300 + 100 * i if c.get("temperatures") else 300
            eps = (
                (0.1 * (i + 1))
                if c.get("hrex") and not c.get("sameHamiltonians")
                else 0.1
            )
            sep = 0.45 + 0.1 * i
            top = f"""[ defaults ]
1 2 no
[ atomtypes ]
Ar 39.948 0.0 A 0.34 {eps}
[ moleculetype ]
AR 1
[ atoms ]
1 Ar 1 AR AR 1 0.0 39.948
[ system ]
Two argon atoms
[ molecules ]
AR 2
"""
            (d / "topol.top").write_text(top)
            (d / "conf.gro").write_text(
                "Argon\n    2\n"
                + f"{1:5d}{'AR':<5}{'AR':>5}{1:5d}{1.0:8.3f}{1.0:8.3f}{1.0:8.3f}\n"
                + f"{2:5d}{'AR':<5}{'AR':>5}{2:5d}{1 + sep:8.3f}{1.0:8.3f}{1.0:8.3f}\n"
                + "3.0 3.0 3.0\n"
            )
            mdp = f"""integrator = md
nsteps = 31
dt = 0.0001
nstlist = 5
cutoff-scheme = Verlet
rlist = 1.05
verlet-buffer-tolerance = -1
coulombtype = Cut-off
rcoulomb = 1.0
vdwtype = Cut-off
vdw-modifier = None
rvdw = 1.0
pbc = xyz
tcoupl = v-rescale
tc-grps = System
tau-t = 0.1
ref-t = {t}
pcoupl = no
constraints = none
gen-vel = yes
gen-temp = {t}
gen-seed = {12345 + i}
ld-seed = {23456 + i}
nstenergy = 1
nstlog = 1
nstcalcenergy = 1
nstcomm = 1
nstxout = 1
nstvout = 1
nstfout = 1
"""
            (d / "md.mdp").write_text(mdp)
            kind = c["bias"]
            at = (
                0.3 + 0.2 * (c["n"] - 1 - i if c.get("reverseBias") else i)
                if c.get("different", True)
                else 0.5
            )
            inp = "CHECKPOINT FILE=plumed.cpt\nd: DISTANCE ATOMS=1,2\n"
            fields = ["d"]
            if kind in ["energy", "combined", "metad-energy"]:
                inp += "e: ENERGY\n"
                fields += ["e"]
            if kind in ["coordinate", "combined"]:
                inp += f"r: RESTRAINT ARG=d AT={at} KAPPA={c.get('kappa', 100)}\n"
                fields += ["r.bias"]
            if kind in ["energy", "combined", "metad-energy"]:
                inp += f"er: RESTRAINT ARG=e AT={0.05 * i} KAPPA=100\n"
                fields += ["er.bias"]
            if kind in ["metad", "metad-energy"]:
                argument = "e" if kind == "metad-energy" else "d"
                inp += f"m: METAD ARG={argument} SIGMA=0.1 HEIGHT={0.2 * (i + 1)} PACE={c.get('pace', 5)} FILE=HILLS\n"
                fields += ["m.bias"]
            if c.get("custom"):
                inp += f"RANDOM_EXCHANGES SEED={7641 + i if c.get('seedMismatch') else 7641}\n"
            inp += "PRINT ARG=" + ",".join(fields) + " FILE=COLVAR STRIDE=1 FMT=%.12f\n"
            (d / "plumed.dat").write_text(inp)
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
                d,
                "grompp.log",
            )
        return root

    def physical(self, c, i, x):
        eps = 0.1 * (i + 1) if c.get("hrex") and not c.get("sameHamiltonians") else 0.1
        z = (0.34 / x) ** 6
        return 4 * eps * (z * z - z)

    def bias(self, c, i, x, step, hills, post=None):
        if post is None:
            post = getattr(self.args, "bias_order", "pre") == "post"
        kind = c["bias"]
        v = 0.0
        u = self.physical(c, i, x)
        if kind in ["coordinate", "combined"]:
            at = (
                0.3 + 0.2 * (c["n"] - 1 - i if c.get("reverseBias") else i)
                if c.get("different", True)
                else 0.5
            )
            v += 0.5 * c.get("kappa", 100) * (x - at) ** 2
        if kind in ["energy", "combined", "metad-energy"]:
            v += 50 * (u - 0.05 * i) ** 2
        if kind in ["metad", "metad-energy"]:
            for h in hills[i]:
                hs = round(h["time"] / 0.0001)
                if hs < step or (post and hs == step):
                    arg = "e" if kind == "metad-energy" else "d"
                    value = u if kind == "metad-energy" else x
                    v += h["height"] * max(
                        0.0,
                        (
                            math.exp(-0.5 * ((value - h[arg]) / h["sigma_" + arg]) ** 2)
                            - math.exp(-6.25)
                        )
                        / (1 - math.exp(-6.25)),
                    )
        return v

    def verify(self, c, root, prefix="run"):
        rows = [
            table(d)
            for i in range(c["n"])
            for d in [next((root / f"rep{i}").glob("COLVAR*"))]
        ]
        hills = [
            table(next((root / f"rep{i}").glob("HILLS*")))
            if c["bias"] in ["metad", "metad-energy"]
            else []
            for i in range(c["n"])
        ]
        bystep = [{round(r["time"] / 0.0001): r for r in rr} for rr in rows]
        events = {}
        for i in range(c["n"]):
            log = (root / f"rep{i}/{prefix}.log").read_text()
            for m in re.finditer(
                r"PLUMED_GREX step=(\d+) a=(\d+) b=(\d+) correction=(\S+) delta=(\S+)",
                log,
            ):
                step, a, b = int(m[1]), int(m[2]), int(m[3])
                events[(step, a, b)] = (float(m[4]), float(m[5]))
        if not events:
            raise AssertionError(
                "no high-precision GREX diagnostics; use the patched build"
            )
        checks = []
        for (step, a, b), (corr, total) in sorted(events.items()):
            xa, xb = bystep[a][step]["d"], bystep[b][step]["d"]
            ba = 1 / (R * (300 + 100 * a if c.get("temperatures") else 300))
            bb = 1 / (R * (300 + 100 * b if c.get("temperatures") else 300))
            vb = ba * (
                self.bias(c, a, xb, step, hills) - self.bias(c, a, xa, step, hills)
            ) + bb * (
                self.bias(c, b, xa, step, hills) - self.bias(c, b, xb, step, hills)
            )
            ua, ub = self.physical(c, a, xa), self.physical(c, b, xb)
            phys = ba * (self.physical(c, a, xb) - ua) + bb * (
                self.physical(c, b, xa) - ub
            )
            expectedcorr = vb + phys if c.get("hrex") else vb
            expected = vb + phys
            err = max(abs(expectedcorr - corr), abs(expected - total))
            tol = self.args.atol + self.args.rtol * max(
                abs(expected), abs(expectedcorr)
            )
            if err > tol:
                raise AssertionError(
                    f"step{step} pair{a},{b}: correction={corr:.9g}, expected={expectedcorr:.9g}; delta={total:.9g}, expected={expected:.9g}; error={err:.3g} tolerance={tol:.3g}"
                )
            if c["bias"] in ["energy", "combined", "metad-energy"]:
                for i, x in [(a, xa), (b, xb)]:
                    if abs(bystep[i][step]["e"] - self.physical(c, i, x)) > tol:
                        raise AssertionError(
                            "ENERGY disagrees with independent LJ energy"
                        )
            segment = (
                (root / f"rep{a}/{prefix}.log")
                .read_text()
                .split(f"Replica exchange at step {step} ")[-1]
                .split("Replica exchange at step ")[0]
            )
            accepted = bool(
                re.search(
                    r"(?<!\d)" + str(a) + r"\s+x\s+" + str(b) + r"(?!\d)", segment
                )
            )
            if total <= 0 and not accepted:
                raise AssertionError("Nonpositive exponent was rejected")
            if total > 100 and accepted:
                raise AssertionError("Exponent beyond cutoff was accepted")
            # In this tiny system a single 0.1-fs MD step changes distance by
            # far less than 0.01 nm. The replicas start 0.1 nm apart. Check the
            # actual moved configuration, not merely a printed acceptance flag.
            for receiver, donor in [
                (a, b if accepted else a),
                (b, a if accepted else b),
            ]:
                if step + 1 in bystep[receiver] and (
                    abs(bystep[receiver][step + 1]["d"] - bystep[donor][step]["d"])
                    > 0.01
                ):
                    raise AssertionError(
                        f"Wrong configuration moved for pair{a},{b} at step{step}"
                    )
            checks.append(
                {
                    "step": step,
                    "a": a,
                    "b": b,
                    "correction": corr,
                    "expected": expected,
                    "delta": total,
                    "error": err,
                    "accepted": accepted,
                    "probability": min(1, math.exp(-min(700, total))),
                }
            )
        return checks

    def case(self, c):
        start = time.time()
        item = {"case": c}
        try:
            root = self.setup(c)
            n = c["n"] * c["ranks"]
            cmd = (
                self.launch(n)
                + ["-multidir"]
                + [f"rep{i}" for i in range(c["n"])]
                + [
                    "-s",
                    "topol.tpr",
                    "-deffnm",
                    "run",
                    "-plumed",
                    "plumed.dat",
                    "-replex",
                    "5",
                    "-reseed",
                    "4321",
                    "-pin",
                    "off",
                    "-notunepme",
                ]
            )
            if c.get("no_plumed"):
                idx = cmd.index("-plumed")
                del cmd[idx : idx + 2]
            if c.get("hrex"):
                cmd += ["-hrex"]
            if c.get("tuning"):
                cmd.remove("-notunepme")
                cmd += ["-tunepme"]
            if c.get("nex"):
                cmd += ["-nex", "2"]
            extra = shlex.split(self.args.mdrun_args)
            cmd += extra
            if (
                "-pme" in extra
                and extra[extra.index("-pme") + 1] == "gpu"
                and "-npme" not in extra
            ):
                cmd += ["-npme", "1" if c["ranks"] > 1 else "0"]
            if c.get("restart"):
                # The first invocation stops before a full step budget, then continues the
                # same TPR with a checkpoint on an exchange boundary (step10).
                self.run(
                    cmd + ["-nsteps", str(c.get("restartAt", 10))], root, "first.log"
                )
                self.run(cmd + ["-cpi", "run.cpt", "-append"], root, "continued.log")
            else:
                code = self.run(
                    cmd,
                    root,
                    "mdrun.log",
                    check=not (c.get("nex") or c.get("expected_failure")),
                )
                if c.get("expected_failure"):
                    diagnostics = failure_diagnostics(root)
                    if not code or c["expected_failure"] not in diagnostics:
                        raise AssertionError("expected diagnostic not observed")
                    item["checks"] = "expected failure verified"
                    item["status"] = "pass"
                    return item
                if c.get("nex"):
                    if (
                        not code
                        or "multiple exchanges" not in (root / "mdrun.log").read_text()
                    ):
                        raise AssertionError("expected -nex rejection")
                    item["checks"] = "expected failure verified"
                    item["status"] = "pass"
                    return item
            # Foreign evaluations must not produce extra PRINT rows or deposit hills.
            for i in range(c["n"]):
                rows = table(next((root / f"rep{i}").glob("COLVAR*")))
                if [round(row["time"] / 0.0001) for row in rows] != list(range(32)):
                    raise AssertionError("missing or duplicated ordinary PLUMED output")
                if c["bias"] in ["metad", "metad-energy"]:
                    hills = table(next((root / f"rep{i}").glob("HILLS*")))
                    if len(hills) != 31 // c.get("pace", 5):
                        raise AssertionError("unexpected deposition count")
                    bystep = {round(row["time"] / 0.0001): row for row in rows}
                    argument = "e" if c["bias"] == "metad-energy" else "d"
                    for hill in hills:
                        step = round(hill["time"] / 0.0001)
                        if abs(hill[argument] - bystep[step][argument]) > 2e-7:
                            raise AssertionError(
                                "Hill deposited at a foreign or stale CV value"
                            )
            item["checks"] = self.verify(c, root)
            if c.get("restart"):
                reference_case = dict(c, name=c["name"] + "-continuous", restart=False)
                reference = self.setup(reference_case)
                self.run(cmd, reference, "mdrun.log")
                self.verify(c, reference)
                maxerr = 0.0
                for i in range(c["n"]):
                    x = table(next((root / f"rep{i}").glob("COLVAR*")))
                    y = table(next((reference / f"rep{i}").glob("COLVAR*")))
                    if len(x) != len(y):
                        raise AssertionError("restart lost/duplicated output rows")
                    for a, b in zip(x, y):
                        for k in a:
                            error = abs(a[k] - b[k])
                            maxerr = max(maxerr, error)
                            if error > self.args.atol + self.args.rtol * abs(b[k]):
                                raise AssertionError(
                                    f"restart differs from continuous run: {k} {a[k]} vs {b[k]}"
                                )
                item["restart_max_error"] = maxerr
            if c.get("compareGrex"):
                reference_case = dict(
                    c, name=c["name"] + "-grex", hrex=False, compareGrex=False
                )
                reference = self.setup(reference_case)
                reference_cmd = list(cmd)
                reference_cmd.remove("-hrex")
                self.run(reference_cmd, reference, "mdrun.log")
                self.verify(reference_case, reference)
                maxima = {"x": 0.0, "v": 0.0, "f": 0.0}
                budgets = {"x": 0.0, "v": 0.0, "f": 0.0}
                for i in range(c["n"]):
                    extracted = []
                    for directory in [root, reference]:
                        d = directory / f"rep{i}"
                        self.run(
                            [self.args.gmx, "dump", "-f", "run.trr"],
                            d,
                            "trajectory.dump",
                            separate_stderr=True,
                        )
                        values = {k: [] for k in maxima}
                        for match in re.finditer(
                            r"\b([xvf])\[\s*\d+\]\s*=\{([^}]+)\}",
                            (d / "trajectory.dump").read_text(),
                        ):
                            values[match[1]].extend(
                                float(z) for z in match[2].split(",")
                            )
                        if not all(values.values()):
                            raise AssertionError(
                                "missing position/velocity/force trajectory data"
                            )
                        extracted.append(values)
                    for k, maximum in maxima.items():
                        a, b = extracted[0][k], extracted[1][k]
                        if len(a) != len(b):
                            raise AssertionError("HREX/GREX trajectory lengths differ")
                        # Use the vector infinity norm, including for force components
                        # near zero after cancellation in the single-precision water system.
                        tolerance = self.args.atol + self.args.rtol * max(
                            map(abs, a + b)
                        )
                        budgets[k] = max(budgets[k], tolerance)
                        for x, y in zip(a, b):
                            error = abs(x - y)
                            maxima[k] = max(maximum, error)
                            if error > tolerance:
                                raise AssertionError(
                                    f"HREX trial altered {k}: {x} vs {y}; tolerance={tolerance}"
                                )
                item["hrex_grex_max_errors"] = maxima
                item["hrex_grex_tolerances"] = budgets
            item["status"] = "pass"
        except (
            AssertionError,
            RuntimeError,
            OSError,
            ValueError,
            KeyError,
            IndexError,
            subprocess.SubprocessError,
        ) as e:
            item["status"] = "fail"
            item["error"] = str(e)
        finally:
            item["seconds"] = round(time.time() - start, 3)
        return item

    def main(self):
        cases = []
        for n in [2, 3, 4]:
            for ranks in [1, 2]:
                for temps in [False, True]:
                    for bias in [
                        "coordinate",
                        "energy",
                        "combined",
                        "metad",
                        "metad-energy",
                    ]:
                        for hrex in [False, True] if self.args.hrex else [False]:
                            cases.append(
                                {
                                    "name": f"n{n}-r{ranks}-t{int(temps)}-{bias}-h{int(hrex)}",
                                    "n": n,
                                    "ranks": ranks,
                                    "temperatures": temps,
                                    "bias": bias,
                                    "hrex": hrex,
                                }
                            )
        for bias in ["coordinate", "energy", "metad", "metad-energy"]:
            for at in [5, 7, 10]:
                for h in [False, True] if self.args.hrex else [False]:
                    cases.append(
                        {
                            "name": f"restart-{bias}-at{at}-h{int(h)}",
                            "n": 2,
                            "ranks": 1,
                            "bias": bias,
                            "restart": True,
                            "restartAt": at,
                            "hrex": h,
                        }
                    )
        cases += [
            {"name": "zero-bias", "n": 2, "ranks": 1, "bias": "none"},
            {
                "name": "large-positive",
                "n": 2,
                "ranks": 1,
                "bias": "coordinate",
                "kappa": 50000,
            },
            {
                "name": "large-negative",
                "n": 2,
                "ranks": 1,
                "bias": "coordinate",
                "kappa": 50000,
                "reverseBias": True,
            },
            {
                "name": "metad-noncoincident",
                "n": 3,
                "ranks": 2,
                "bias": "metad",
                "pace": 7,
            },
            {
                "name": "identical-bias",
                "n": 2,
                "ranks": 1,
                "bias": "coordinate",
                "different": False,
            },
            {
                "name": "invalid-nex",
                "n": 2,
                "ranks": 1,
                "bias": "coordinate",
                "nex": True,
            },
        ]
        if self.args.custom:
            for h in [False, True] if self.args.hrex else [False]:
                cases.append(
                    {
                        "name": f"custom-h{int(h)}",
                        "n": 4,
                        "ranks": 1,
                        "bias": "coordinate",
                        "custom": True,
                        "hrex": h,
                    }
                )
        if self.args.custom:
            for h in [False, True] if self.args.hrex else [False]:
                for at in [5, 10]:
                    cases.append(
                        {
                            "name": f"restart-custom-at{at}-h{int(h)}",
                            "n": 4,
                            "ranks": 1,
                            "bias": "coordinate",
                            "custom": True,
                            "hrex": h,
                            "restart": True,
                            "restartAt": at,
                        }
                    )
        if self.args.custom:
            cases.append(
                {
                    "name": "invalid-pattern-mismatch",
                    "n": 4,
                    "ranks": 1,
                    "bias": "coordinate",
                    "custom": True,
                    "seedMismatch": True,
                    "expected_failure": "PLUMED exchange permutation",
                }
            )
        if self.args.hrex:
            cases += [
                {
                    "name": "invalid-hrex-without-plumed",
                    "n": 2,
                    "ranks": 1,
                    "bias": "coordinate",
                    "hrex": True,
                    "no_plumed": True,
                    "expected_failure": "requires -plumed",
                },
                {
                    "name": "invalid-hrex-autotuning",
                    "n": 2,
                    "ranks": 1,
                    "bias": "coordinate",
                    "hrex": True,
                    "tuning": True,
                    "expected_failure": "requires -notunepme",
                },
            ]
        if self.args.hrex:
            for n in [2, 3, 4]:
                for ranks in [1, 2]:
                    for temps in [False, True]:
                        cases.append(
                            {
                                "name": f"physical-only-n{n}-r{ranks}-t{int(temps)}",
                                "n": n,
                                "ranks": ranks,
                                "temperatures": temps,
                                "bias": "none",
                                "hrex": True,
                            }
                        )
        if self.args.custom:
            for bias in ["combined", "metad", "metad-energy"]:
                for h in [False, True] if self.args.hrex else [False]:
                    for at in [7, 10]:
                        cases.append(
                            {
                                "name": f"custom-stress-{bias}-at{at}-h{int(h)}",
                                "n": 3,
                                "ranks": 2,
                                "temperatures": True,
                                "bias": bias,
                                "custom": True,
                                "hrex": h,
                                "restart": True,
                                "restartAt": at,
                            }
                        )
        if self.args.hrex:
            for bias in ["coordinate", "energy", "combined", "metad", "metad-energy"]:
                cases.append(
                    {
                        "name": f"hrex-grex-equivalence-{bias}",
                        "n": 3,
                        "ranks": 2,
                        "temperatures": True,
                        "bias": bias,
                        "hrex": True,
                        "custom": self.args.custom,
                        "sameHamiltonians": True,
                        "compareGrex": True,
                    }
                )
        for c in cases:
            mode = getattr(self.args, "exchange_mode", "both")
            if mode == "bias" and c.get("hrex"):
                continue
            if mode == "hrex" and not c.get("hrex"):
                continue
            if self.args.filter and self.args.filter not in c["name"]:
                continue
            item = self.case(c)
            self.results.append(item)
            print(json.dumps(item), flush=True)
            (self.root / "results.json").write_text(
                json.dumps(self.results, indent=2) + "\n"
            )
        if not self.results:
            raise SystemExit("No test cases selected")
        return int(any(x["status"] == "fail" for x in self.results))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gmx", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--serial-launcher", default="")
    p.add_argument("--launcher", default="mpiexec -n {n}")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--mdrun-args", default="")
    p.add_argument("--timeout", type=int, default=120)
    p.add_argument("--atol", type=float, default=3e-5)
    p.add_argument("--rtol", type=float, default=3e-5)
    p.add_argument("--exchange-mode", choices=["both", "bias", "hrex"], default="both")
    p.add_argument("--bias-order", choices=["pre", "post"], default="pre")
    p.add_argument("--filter", default="")
    p.add_argument("--hrex", action="store_true")
    p.add_argument("--custom", action="store_true")
    a = p.parse_args()
    a.gmx = str(Path(a.gmx).resolve())
    sys.exit(Suite(a).main())
