#!/usr/bin/env python3
"""Check the installed multi-target OpenACC library on host or NVIDIA hardware."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path, required=True)
parser.add_argument("--device", choices=("host", "nvidia"), default="host")
args = parser.parse_args()
prefix = args.prefix.resolve()
source = Path(__file__).resolve().parents[2]
program = str(prefix / "bin/plumed")
environment = dict(
    os.environ,
    PATH=str(prefix / "bin") + os.pathsep + os.environ.get("PATH", ""),
    PLUMED_AUTOLOAD="1",
    PLUMED_NUM_THREADS="2",
    OMP_NUM_THREADS="2",
    PLUMED_KERNEL=str(prefix / "lib/libplumedKernel.so"),
)
environment["NVCOMPILER_ACC_DEVICE_TYPE"] = args.device.upper()
environment["ACC_DEVICE_TYPE"] = args.device.upper()
if args.device == "nvidia":
    subprocess.run(["nvidia-smi"], check=True, timeout=30)
    environment["NVCOMPILER_ACC_NOTIFY"] = "1"
for mode in ("autoload", "explicit"):
    env = dict(environment, PLUMED_AUTOLOAD="1" if mode == "autoload" else "0")
    command = [program] + ([] if mode == "autoload" else ["--load-module=openacc_ptm"])
    subprocess.run(
        [*command, "modules", "-q", "--has", "openacc_ptm"],
        env=env,
        check=True,
        timeout=30,
    )
    subprocess.run(
        [*command, "modules", "-q", "--has-feature", "openacc"],
        env=env,
        check=True,
        timeout=30,
    )
    for action in ("DISTANCE",):
        result = subprocess.run(
            [*command, "manual", "--action", action],
            env=env,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert "USEGPU" in result.stdout, action
assert (
    subprocess.run(
        [program, "modules", "-q", "--has", "openacc_ptm"],
        env=dict(environment, PLUMED_AUTOLOAD="0"),
        check=False,
        timeout=30,
    ).returncode
    == 1
)
with tempfile.TemporaryDirectory(prefix="plumed-openacc-regressions-") as scratch:
    regtest = Path(scratch) / "regtest"
    regtest.mkdir()
    (regtest / "scripts").symlink_to(
        source / "regtest/scripts", target_is_directory=True
    )
    (regtest / "trajectories").symlink_to(
        source / "plugins/openaccPTM/regtest/trajectories", target_is_directory=True
    )
    for mode in ("autoload", "explicit"):
        for group, name in (
            ("targetGPU", "rt-GPU-DISTANCE"),
            ("targetSAXS", "rt-GPU-saxs-chi2-deriv"),
        ):
            fixture = regtest / mode / name
            shutil.copytree(
                source / "plugins/openaccPTM/regtest" / group / name, fixture
            )
            if group == "targetSAXS":
                # Relocate the original fixture's shared CPU references/input.
                for path in (source / "regtest/isdb/rt-saxs-chi2-deriv").iterdir():
                    if path.name.endswith(".reference") or path.name == "twoatom.pdb":
                        shutil.copy2(path, fixture / path.name)
                config = fixture / "config"
                config.write_text(
                    "\n".join(
                        line
                        for line in config.read_text().splitlines()
                        if not line.startswith("extra_files=")
                    )
                    + "\n"
                )
            input_path = fixture / "plumed.dat"
            input_path.write_text(
                input_path.read_text().replace(
                    "LOAD FILE=../../../../plumedOpenACC.so",
                    "LOAD MODULE=openacc_ptm" if mode == "explicit" else "",
                )
            )
            result = subprocess.run(
                ["bash", str(source / "regtest/scripts/run")],
                cwd=fixture,
                env=dict(
                    environment,
                    PLUMED_AUTOLOAD="0" if mode == "explicit" else "1",
                    PLUMED_PROGRAM_NAME=program,
                ),
                capture_output=True,
                text=True,
                check=False,
                timeout=120,
            )
            driver_log = "\n".join(
                path.read_text(errors="replace")
                for path in (fixture / "tmp/out", fixture / "tmp/err")
                if path.exists()
            )
            report = (
                (fixture / "report.txt").read_text()
                if (fixture / "report.txt").exists()
                else ""
            )
            if (
                result.returncode
                or not report
                or any(
                    marker in report
                    for marker in (
                        "FAILURE",
                        "NOT_APPLICABLE",
                        "DEFERRED",
                        "SKIP_ON_TRAVIS",
                    )
                )
            ):
                raise AssertionError(
                    f"{mode} {name}:\n{report}\n{result.stdout}\n{result.stderr}\n{driver_log}"
                )
            if args.device == "nvidia":
                assert "launch CUDA kernel" in driver_log, (
                    "No GPU kernel launch observed"
                )
            print(
                f"Installed OpenACC numerical regression passed: {args.device} {mode} {name}",
                flush=True,
            )
