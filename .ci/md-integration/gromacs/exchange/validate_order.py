"""Count real forward calculations, force calls and check local force/deposition state."""

import argparse
from collections import Counter
import importlib.util
import json
import itertools
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from types import SimpleNamespace

location = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("analytic", location / "validate.py")
analytic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analytic)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--kernel", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--launcher", required=True)
    options = parser.parse_args()
    root = Path(options.output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    os.environ["PLUMED_KERNEL"] = str(Path(options.kernel).resolve())
    prefix = Path(options.kernel).resolve().parent.parent
    os.environ["PATH"] = str(prefix / "bin") + os.pathsep + os.environ["PATH"]
    shutil.copy2(location / "count-exchange.cpp", root / "count-exchange.cpp")
    with (root / "mklib.log").open("w") as log:
        subprocess.run(
            [str(prefix / "bin/plumed"), "mklib", "count-exchange.cpp"],
            cwd=root,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    library = root / (
        "count-exchange.dylib" if sys.platform == "darwin" else "count-exchange.so"
    )
    args = SimpleNamespace(
        gmx=str(Path(options.gmx).resolve()),
        output=str(root),
        threads=1,
        timeout=120,
        serial_launcher="",
        launcher=options.launcher,
        mdrun_args="-nb cpu -pme cpu -update cpu",
        atol=3e-5,
        rtol=3e-5,
    )

    class Suite(analytic.Suite):
        def setup(self, case):
            result = super().setup(case)
            for replica in range(case["n"]):
                path = result / f"rep{replica}/plumed.dat"
                text = path.read_text()
                text = f"LOAD FILE={library}\n" + text
                argument = "e" if "energy" in case["bias"] else "d"
                line = "e: ENERGY\n" if argument == "e" else "d: DISTANCE ATOMS=1,2\n"
                text = text.replace(
                    line, line + f"counted: COUNT_EXCHANGE ARG={argument}\n"
                )
                text = text.replace(
                    f"r: RESTRAINT ARG={argument}", "r: RESTRAINT ARG=counted"
                )
                text = text.replace(f"METAD ARG={argument}", "METAD ARG=counted")
                text = re.sub(
                    r"((?:r: RESTRAINT|m: METAD)[^\n]+)",
                    rf"\1 STRIDE={case['stride']}",
                    text,
                )
                text = re.sub(r"PRINT ARG=\S+", "PRINT ARG=d", text)
                path.write_text(text)
            return result

    suite = Suite(args)
    results = []
    for hrex in [False, True]:
        for kind, stride in itertools.product(
            ["coordinate", "energy", "metad", "metad-energy"], [1, 4, 10]
        ):
            case = dict(
                name=f"count-{kind}-h{int(hrex)}-stride{stride}",
                n=2,
                ranks=1,
                bias=kind,
                hrex=hrex,
                stride=stride,
                pace=5 if stride == 1 else 2 * stride,
            )
            directory = suite.setup(case)
            command = suite.launch(2) + [
                "-multidir",
                "rep0",
                "rep1",
                "-s",
                "topol.tpr",
                "-deffnm",
                "run",
                "-plumed",
                "plumed.dat",
                "-replex",
                "5",
                "-reseed",
                "2026",
                "-notunepme",
                "-pin",
                "off",
            ]
            command += shlex.split(args.mdrun_args)
            if hrex:
                command += ["-hrex"]
            suite.run(command, directory, "mdrun.log")
            for replica in range(2):
                local = directory / f"rep{replica}"
                log = next(local.glob("PLUMED.OUT*")).read_text()
                counts = Counter(
                    (int(step), int(exchange))
                    for step, exchange in re.findall(
                        r"EXCHANGE_COUNT step=(\d+) exchange=\d+ auxiliary=(\d+)", log
                    )
                )
                expected = Counter()
                for step in range(32):
                    if step % stride == 0:
                        expected[step, 0] += 1
                    if step in range(5, 31, 5):
                        expected[step, 1] += 2
                if counts != expected:
                    raise AssertionError(
                        f"{case['name']}: forward calls {counts} expected {expected}"
                    )
                force_count = re.search(
                    r"^ Force\s+\d+\s+\d+\s+(\d+)\s",
                    (local / "run.log").read_text(),
                    re.M,
                )
                if force_count is None or int(force_count[1]) != 32 + 6 * hrex:
                    raise AssertionError("Unexpected physical force evaluation count")
                # The full force trajectory must retain normal MTS scaling.
                suite.run(
                    [args.gmx, "dump", "-f", "run.trr"],
                    local,
                    "trajectory.dump",
                    separate_stderr=True,
                )
                frames = re.split(
                    r"(?=^.*?frame\s+\d+:)",
                    (local / "trajectory.dump").read_text(),
                    flags=re.M,
                )
                distances = {
                    round(r["time"] / 0.0001): r["d"]
                    for r in analytic.table(next(local.glob("COLVAR*")))
                }
                hills = (
                    analytic.table(next(local.glob("HILLS*")))
                    if "metad" in case["bias"]
                    else []
                )
                if "metad" in case["bias"]:
                    if len(hills) != 31 // case["pace"]:
                        raise AssertionError("Unexpected hill count with MTS stride")
                    for hill in hills:
                        hill_step = round(hill["time"] / 0.0001)
                        center = distances[hill_step]
                        if "energy" in case["bias"]:
                            z = (0.34 / center) ** 6
                            center = (
                                4 * (0.1 * (replica + 1) if hrex else 0.1) * (z * z - z)
                            )
                        if abs(hill["counted"] - center) > 3e-7:
                            raise AssertionError(
                                "MTS hill deposited at foreign/stale coordinates"
                            )
                checked = 0
                for frame in frames:
                    step_match = re.search(r"\bstep=\s*(\d+)", frame)
                    if not step_match:
                        continue
                    step = int(step_match[1])
                    vectors = {}
                    for kind, atom, values in re.findall(
                        r"\b([xf])\[\s*(\d+)\]\s*=\{([^}]+)\}", frame
                    ):
                        vectors[kind, int(atom)] = list(map(float, values.split(",")))
                    if len(vectors) != 4:
                        raise AssertionError("Missing trajectory positions/forces")
                    dx = [b - a for a, b in zip(vectors["x", 0], vectors["x", 1])]
                    direction_norm = math.sqrt(sum(x * x for x in dx))
                    distance = distances[step]
                    z = (0.34 / distance) ** 6
                    eps = 0.1 * (replica + 1) if hrex else 0.1
                    derivative = 24 * eps * (z - 2 * z * z) / distance
                    if step % stride == 0:
                        energy = 4 * eps * (z * z - z)
                        argument = energy if "energy" in case["bias"] else distance
                        bias_derivative = 0.0
                        if case["bias"] == "coordinate":
                            bias_derivative += 100 * (distance - (0.3 + 0.2 * replica))
                        if "energy" in case["bias"]:
                            bias_derivative += 100 * (energy - 0.05 * replica)
                        for hill in hills:
                            if round(hill["time"] / 0.0001) >= step:
                                continue
                            difference = argument - hill["counted"]
                            sigma = hill["sigma_counted"]
                            scaled = difference / sigma
                            if scaled * scaled < 12.5:
                                bias_derivative -= (
                                    hill["height"]
                                    * difference
                                    / (sigma * sigma)
                                    * math.exp(-0.5 * scaled * scaled)
                                    / (1 - math.exp(-6.25))
                                )
                        if "energy" in case["bias"]:
                            derivative *= 1 + stride * bias_derivative
                        else:
                            derivative += stride * bias_derivative
                    for atom in [0, 1]:
                        expected_force = [
                            (1 if atom == 0 else -1) * derivative * x / direction_norm
                            for x in dx
                        ]
                        error = max(
                            abs(a - b)
                            for a, b in zip(vectors["f", atom], expected_force)
                        )
                        if error > 0.003 + 1e-4 * max(map(abs, expected_force)):
                            raise AssertionError(
                                f"Wrong MTS force at step {step}: {vectors['f', atom]} expected {expected_force}"
                            )
                    checked += 1
                if checked != 32:
                    raise AssertionError(
                        f"Checked {checked} force frames instead of 32"
                    )
            item = dict(
                case=case, status="pass", force_frames=64, forward_counts="exact"
            )
            results.append(item)
            (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
            print(case["name"], "pass", flush=True)
    # Native lambda cross energies must not trigger another physical force call.
    lambda_spec = importlib.util.spec_from_file_location(
        "lambda_suite", location / "validate_lambda.py"
    )
    lambda_module = importlib.util.module_from_spec(lambda_spec)
    lambda_spec.loader.exec_module(lambda_module)
    lambda_suite = lambda_module.Suite(args)
    for bias in ["coordinate", "energy", "metad-energy"]:
        case = dict(name=f"native-lambda-force-count-{bias}", n=2, ranks=1, bias=bias)
        item = lambda_suite.case(case)
        if item["status"] != "pass":
            raise AssertionError(item)
        for replica in range(2):
            text = (root / case["name"] / f"rep{replica}/run.log").read_text()
            count = re.search(r"^ Force\s+\d+\s+\d+\s+(\d+)\s", text, re.M)
            if count is None or int(count[1]) != 32:
                raise AssertionError(
                    "Native lambda exchange performed an extra force-field calculation"
                )
        item["force_calls_per_replica"] = 32
        results.append(item)
        (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
        print(case["name"], "pass", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
