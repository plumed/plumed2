#!/usr/bin/env python3
"""Check installed CUDA module loading; GPU numerical tests are explicit."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path, required=True)
parser.add_argument("--cuda-version", required=True)
parser.add_argument(
    "--gpu", action="store_true", help="require and execute GPU numerical regressions"
)
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
)
registration = (prefix / "lib/plumed/plugins/cuda_coordination.plumed").read_text()
assert "variant cuda" + args.cuda_version + "\n" in registration
for mode in ("autoload", "explicit"):
    env = dict(environment, PLUMED_AUTOLOAD="1" if mode == "autoload" else "0")
    command = [program] + (
        [] if mode == "autoload" else ["--load-module=cuda_coordination"]
    )
    subprocess.run(
        [*command, "modules", "-q", "--has", "cuda_coordination"],
        env=env,
        check=True,
        timeout=30,
    )
    subprocess.run(
        [*command, "modules", "-q", "--has-feature", "cuda"],
        env=env,
        check=True,
        timeout=30,
    )
    for action in ("CUDACOORDINATION", "CUDACOORDINATIONFLOAT"):
        result = subprocess.run(
            [*command, "manual", "--action", action],
            env=env,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert "GROUPA" in result.stdout and "R_0" in result.stdout, action
assert (
    subprocess.run(
        [program, "modules", "-q", "--has", "cuda_coordination"],
        env=dict(environment, PLUMED_AUTOLOAD="0"),
        check=False,
        timeout=30,
    ).returncode
    == 1
)
print(
    f"CUDA {args.cuda_version}: native loading, actions and compiled feature checks passed",
    flush=True,
)
if not args.gpu:
    print(
        "GPU numerical tests not requested; no GPU calculation has been validated",
        flush=True,
    )
    raise SystemExit(0)

with tempfile.TemporaryDirectory(prefix="plumed-cuda-regressions-") as scratch:
    regtest = Path(scratch) / "regtest"
    regtest.mkdir()
    (regtest / "scripts").symlink_to(
        source / "regtest/scripts", target_is_directory=True
    )
    (regtest / "trajectories").symlink_to(
        source / "plugins/cudaCoord/regtest/trajectories", target_is_directory=True
    )
    for mode in ("autoload", "explicit"):
        for name in ("rt-double-nopbc", "rt-float-nopbc"):
            fixture = regtest / mode / name
            shutil.copytree(
                source / "plugins/cudaCoord/regtest/cudatest" / name, fixture
            )
            input_path = fixture / "plumed.dat"
            input_path.write_text(
                input_path.read_text().replace(
                    "LOAD FILE=../../../../CudaCoordination.so",
                    "LOAD MODULE=cuda_coordination" if mode == "explicit" else "",
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
            print(
                f"Installed CUDA numerical regression passed: {mode} {name}", flush=True
            )
