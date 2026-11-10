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
        for name in ("rt-GPU-DISTANCE",):
            fixture = regtest / mode / name
            shutil.copytree(
                source / "plugins/openaccPTM/regtest/targetGPU" / name, fixture
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
                    f"{mode} {name}:\n{report}\n{result.stdout}\n{result.stderr}"
                )
            if args.device == "nvidia":
                log = (
                    (fixture / "tmp/err").read_text()
                    if (fixture / "tmp/err").exists()
                    else ""
                )
                assert "launch CUDA kernel" in result.stderr + result.stdout + log, (
                    "No GPU kernel launch observed"
                )
            print(
                f"Installed OpenACC numerical regression passed: {args.device} {mode} {name}",
                flush=True,
            )
