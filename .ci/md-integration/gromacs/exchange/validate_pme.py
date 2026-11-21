#!/usr/bin/env python3
"""PME oracle: independent CPU reruns under every Hamiltonian on every trajectory.
The measured acceptance exponent must match cross energies; trajectories and
random decisions need not agree across CPU/GPU architectures.
"""

import importlib.util
import sys
import json
import re
import shlex
import math
import argparse
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "validate", Path(__file__).with_name("validate.py")
)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


class PMESuite(v.Suite):
    def setup(self, c):
        root = super().setup(c)
        for i in range(c["n"]):
            d = root / f"rep{i}"
            eps = (
                0.1 * (i + 1)
                if c.get("hrex") and not c.get("sameHamiltonians")
                else 0.1
            )
            charge = (
                0.04 * math.sqrt(i + 1)
                if c.get("hrex") and not c.get("sameHamiltonians")
                else 0.04
            )
            (d / "topol.top").write_text(f"""[ defaults ]
1 2 no
[ atomtypes ]
Ar 39.948 0.0 A 0.34 {eps}
[ moleculetype ]
AR 0
[ atoms ]
1 Ar 1 AR AR 1 {charge} 39.948
2 Ar 2 AR AR 2 {-charge} 39.948
[ system ]
Neutral pair PME
[ molecules ]
AR 1
""")
            mdp = (
                (d / "md.mdp")
                .read_text()
                .replace(
                    "coulombtype = Cut-off",
                    "coulombtype = PME\nfourier-nx = 32\nfourier-ny = 32\nfourier-nz = 32\npme-order = 4\newald-rtol = 1e-6",
                )
            )
            if c.get("mismatch") and i == 1:
                mdp = mdp.replace("fourier-nx = 32", "fourier-nx = 36")
            if c.get("water"):
                (d / "conf.gro").write_text(
                    Path(__file__).with_name("spc216.gro").read_text()
                )
                scale = (
                    (1.0 + 0.1 * i)
                    if c.get("hrex") and not c.get("sameHamiltonians")
                    else 1.0
                )
                (d / "topol.top").write_text(f"""[ defaults ]
1 2 no
[ atomtypes ]
OW 15.9994 0.0 A 0.316557 {0.650194 * scale}
HW 1.008 0.0 A 0 0
[ moleculetype ]
SOL 2
[ atoms ]
1 OW 1 SOL OW 1 {-0.8476 * math.sqrt(scale)} 15.9994
2 HW 1 SOL HW1 1 {0.4238 * math.sqrt(scale)} 1.008
3 HW 1 SOL HW2 1 {0.4238 * math.sqrt(scale)} 1.008
[ settles ]
1 1 0.1 0.1633
[ exclusions ]
1 2 3
2 1 3
3 1 2
[ system ]
216 rigid water molecules
[ molecules ]
SOL 216
""")
                mdp = (
                    mdp.replace("rlist = 1.05", "rlist = 0.85")
                    .replace("rcoulomb = 1.0", "rcoulomb = 0.8")
                    .replace("rvdw = 1.0", "rvdw = 0.8")
                )
                inp = (d / "plumed.dat").read_text().replace("ATOMS=1,2", "ATOMS=1,4")
                # Replace specifically the ENERGY restraint even for combined input.
                inp = re.sub(r"(er: RESTRAINT[^\n]*KAPPA=)100", r"\g<1>0.000001", inp)
                (d / "plumed.dat").write_text(inp)
            (d / "md.mdp").write_text(mdp)
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
                "grompp-pme.log",
            )
        return root

    def verify(self, c, root, prefix="run"):
        cross = {}
        layout_error = {}
        noise = []
        repeat_noise = []
        gpu_cpu = []
        for i in range(c["n"]):
            d = root / f"rep{i}"
            for j in range(c["n"]):
                name = f"oracle-{j}"
                cmd = self.launch(1) + [
                    "-s",
                    "topol.tpr",
                    "-rerun",
                    str(root / f"rep{j}/{prefix}.trr"),
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
                self.run(cmd, d, name + ".log")
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
                    ],
                    d,
                    name + "-energy.log",
                    stdin="Potential\n0\n",
                )
                cross[i, j] = {
                    round(float(line.split()[0]) / 0.0001): float(line.split()[1])
                    for line in (d / (name + ".xvg")).read_text().splitlines()
                    if line and not line.startswith(("#", "@"))
                }
                layout_error[i, j] = {step: 0.0 for step in cross[i, j]}
                if c.get("water") or c["ranks"] > 1:
                    rn = f"layout-{j}"
                    cmd_layout = self.launch(c["ranks"]) + [
                        "-s",
                        "topol.tpr",
                        "-rerun",
                        str(root / f"rep{j}/{prefix}.trr"),
                        "-deffnm",
                        rn,
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
                    self.run(cmd_layout, d, rn + ".log")
                    self.run(
                        [
                            self.args.gmx,
                            "energy",
                            "-f",
                            rn + ".edr",
                            "-o",
                            rn + ".xvg",
                            "-xvg",
                            "none",
                        ],
                        d,
                        rn + "-energy.log",
                        stdin="Potential\n0\n",
                    )
                    alternative = {
                        round(float(l.split()[0]) / 0.0001): float(l.split()[1])
                        for l in (d / (rn + ".xvg")).read_text().splitlines()
                        if l and not l.startswith(("#", "@"))
                    }
                    for step in cross[i, j]:
                        spread = abs(alternative[step] - cross[i, j][step])
                        layout_error[i, j][step] = spread
                        if spread > self.args.energy_atol + self.args.rtol * abs(
                            cross[i, j][step]
                        ):
                            raise AssertionError(
                                f"Independent MPI layouts disagree beyond energy budget: {spread}"
                            )
                if getattr(self.args, "gpu_repeats", False):
                    repeats = []
                    for repeat in [1, 2]:
                        rn = f"gpu-repeat-{j}-{repeat}"
                        gpu = (
                            self.launch(1)
                            + [
                                "-s",
                                "topol.tpr",
                                "-rerun",
                                str(root / f"rep{j}/{prefix}.trr"),
                                "-deffnm",
                                rn,
                                "-pin",
                                "off",
                                "-notunepme",
                            ]
                            + shlex.split(self.args.mdrun_args)
                        )
                        self.run(gpu, d, rn + ".log")
                        self.run(
                            [
                                self.args.gmx,
                                "energy",
                                "-f",
                                rn + ".edr",
                                "-o",
                                rn + ".xvg",
                                "-xvg",
                                "none",
                            ],
                            d,
                            rn + "-energy.log",
                            stdin="Potential\n0\n",
                        )
                        repeats.append(
                            {
                                round(float(l.split()[0]) / 0.0001): float(l.split()[1])
                                for l in (d / (rn + ".xvg")).read_text().splitlines()
                                if l and not l.startswith(("#", "@"))
                            }
                        )
                    for step in cross[i, j]:
                        measured = abs(repeats[0][step] - repeats[1][step])
                        repeat_noise.append(measured)
                        difference = max(
                            abs(rr[step] - cross[i, j][step]) for rr in repeats
                        )
                        gpu_cpu.append(difference)
                        layout_error[i, j][step] = max(
                            layout_error[i, j][step], measured, difference
                        )
                        tol = self.args.energy_atol + self.args.rtol * abs(
                            cross[i, j][step]
                        )
                        if max(measured, difference) > tol:
                            raise AssertionError(
                                f"GPU energy noise/reference difference exceeds budget: repeat={measured}, GPU-CPU={difference}, budget={tol}"
                            )

        rows = [
            {
                round(row["time"] / 0.0001): row
                for row in v.table(next((root / f"rep{i}").glob("COLVAR*")))
            }
            for i in range(c["n"])
        ]
        events = {}
        for i in range(c["n"]):
            for m in re.finditer(
                r"PLUMED_GREX step=(\d+) a=(\d+) b=(\d+) correction=(\S+) delta=(\S+)",
                (root / f"rep{i}/{prefix}.log").read_text(),
            ):
                events[tuple(map(int, m.group(1, 2, 3)))] = tuple(
                    map(float, m.group(4, 5))
                )
        if not events:
            raise AssertionError("No GREX events")
        checks = []

        def bias(i, j, step):
            x = rows[j][step]["d"]
            e = cross[i, j][step]
            kind = c["bias"]
            r = 0.0
            if kind in ["coordinate", "combined"]:
                r += 50 * (x - (0.3 + 0.2 * i)) ** 2
            if kind in ["energy", "combined"]:
                r += (0.0000005 if c.get("water") else 50) * (e - 0.05 * i) ** 2
            return r

        for (step, a, b), (corr, delta) in sorted(events.items()):
            ba = 1 / (v.R * (300 + 100 * a if c.get("temperatures") else 300))
            bb = 1 / (v.R * (300 + 100 * b if c.get("temperatures") else 300))
            phys = ba * (cross[a, b][step] - cross[a, a][step]) + bb * (
                cross[b, a][step] - cross[b, b][step]
            )
            bterm = ba * (bias(a, b, step) - bias(a, a, step)) + bb * (
                bias(b, a, step) - bias(b, b, step)
            )
            expected = phys + bterm
            expectedcorr = expected if c.get("hrex") else bterm
            err = max(abs(expected - delta), abs(expectedcorr - corr))
            tol = self.args.atol + self.args.rtol * max(
                abs(expected), abs(expectedcorr)
            )
            # Propagate independently measured MPI reduction variation through
            # physical energy and the ENERGY-bias derivative, not a blanket
            # large-system exponent tolerance. The small analytic cases retain
            # their strict, architecture-independent bounds.
            propagated = 0.0
            ek = 0.000001 if c.get("water") else 100.0
            for i, j, beta in [(a, b, ba), (a, a, ba), (b, a, bb), (b, b, bb)]:
                spread = layout_error[i, j][step]
                derivative = (
                    ek * abs(cross[i, j][step] - 0.05 * i)
                    if c["bias"] in ["energy", "combined"]
                    else 0.0
                )
                propagated += beta * spread * (1.0 + derivative)
            tol += 2 * propagated

            if err > tol:
                raise AssertionError(
                    f"PME rerun exponent mismatch step{step} pair{a},{b}: {delta} expected{expected}, correction{corr} expected{expectedcorr}, error{err} tolerance{tol}"
                )
            for i in [a, b]:
                if "e" in rows[i][step]:
                    diff = abs(rows[i][step]["e"] - cross[i, i][step])
                    noise.append(diff)
                    if diff > self.args.energy_atol + self.args.rtol * abs(
                        cross[i, i][step]
                    ):
                        raise AssertionError(
                            f"Own ENERGY differs from CPU rerun by {diff}"
                        )
            checks.append(
                dict(
                    step=step,
                    a=a,
                    b=b,
                    delta=delta,
                    expected=expected,
                    error=err,
                    tolerance=tol,
                    layout_uncertainty=propagated,
                )
            )
        return dict(
            events=checks,
            max_own_energy_difference=max(noise, default=0.0),
            max_repeat_energy_noise=max(repeat_noise, default=0.0),
            max_gpu_cpu_energy_difference=max(gpu_cpu, default=0.0),
            oracle="independent CPU PME reruns",
        )

    def main(self):
        for n in [2, 3]:
            for ranks in [1, 2]:
                for kind in ["coordinate", "energy", "combined"]:
                    for h in [False, True]:
                        c = dict(
                            name=f"pme-n{n}-r{ranks}-{kind}-h{int(h)}",
                            n=n,
                            ranks=ranks,
                            bias=kind,
                            hrex=h,
                            temperatures=True,
                            custom=self.args.custom,
                        )
                        mode = getattr(self.args, "exchange_mode", "both")
                        if mode == "bias" and h:
                            continue
                        if mode == "hrex" and not h:
                            continue
                        if self.args.filter and self.args.filter not in c["name"]:
                            continue
                        item = self.case(c)
                        self.results.append(item)
                        print(json.dumps(item), flush=True)
                        (self.root / "results.json").write_text(
                            json.dumps(self.results, indent=2) + "\n"
                        )
        if getattr(self.args, "water", False):
            for h in [False, True]:
                for ranks, kind in [
                    (1, "coordinate"),
                    (2, "energy"),
                    (2, "combined"),
                    (2, "coordinate"),
                ]:
                    c = dict(
                        name=f"water-r{ranks}-{kind}-h{int(h)}",
                        n=2,
                        ranks=ranks,
                        bias=kind,
                        hrex=h,
                        temperatures=True,
                        water=True,
                        custom=self.args.custom,
                        sameHamiltonians=(h and ranks == 2 and kind == "coordinate"),
                        compareGrex=(h and ranks == 2 and kind == "coordinate"),
                    )
                    mode = getattr(self.args, "exchange_mode", "both")
                    if mode == "bias" and h:
                        continue
                    if mode == "hrex" and not h:
                        continue
                    if self.args.filter and self.args.filter not in c["name"]:
                        continue
                    item = self.case(c)
                    self.results.append(item)
                    print(json.dumps(item), flush=True)
                    (self.root / "results.json").write_text(
                        json.dumps(self.results, indent=2) + "\n"
                    )
        if (
            not self.args.filter
            and getattr(self.args, "exchange_mode", "both") != "bias"
        ):
            c = dict(
                name="pme-mismatch-rejected",
                n=2,
                ranks=1,
                bias="coordinate",
                hrex=True,
                mismatch=True,
                expected_failure="HREX PME nx",
            )
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
    p.add_argument("--timeout", type=int, default=180)
    p.add_argument("--atol", type=float, default=0.002)
    p.add_argument("--rtol", type=float, default=3e-5)
    p.add_argument("--energy-atol", type=float, default=0.0001)
    p.add_argument("--water", action="store_true")
    p.add_argument("--gpu-repeats", action="store_true")
    p.add_argument("--filter", default="")
    p.add_argument("--custom", action="store_true")
    p.add_argument("--exchange-mode", choices=["both", "bias", "hrex"], default="both")
    a = p.parse_args()
    a.gmx = str(Path(a.gmx).resolve())
    sys.exit(PMESuite(a).main())
