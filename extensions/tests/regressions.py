#!/usr/bin/env python3
"""Run representative existing numerical fixtures against installed packages."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", required=True, type=Path)
args = parser.parse_args()
prefix = args.prefix.resolve()
source = Path(__file__).resolve().parents[2]
fixtures = [
    "ves/rt-le-1d-static-bias",
    "fourier/rt-fftw",
    "drr/rt-checkpoint-drr",
    "pytorch/rt-pytorch_model_strided",
]
environment = dict(
    os.environ,
    PLUMED_PROGRAM_NAME=str(prefix / "bin/plumed"),
    PLUMED_AUTOLOAD="1",
    PLUMED_NUM_THREADS="2",
    OMP_NUM_THREADS="2",
)
with tempfile.TemporaryDirectory(prefix="plumed-installed-regressions-") as scratch:
    regtest = Path(scratch) / "regtest"
    regtest.mkdir()
    for name in ("scripts", "trajectories"):
        (regtest / name).symlink_to(source / "regtest" / name, target_is_directory=True)
    for name in fixtures:
        fixture = regtest / name
        shutil.copytree(
            source / "regtest" / name,
            fixture,
            ignore=shutil.ignore_patterns("tmp", "report.txt"),
        )
        # Original config, feature/module guards, hooks and references stay intact.
        print(f"Running installed numerical regression: {name}", flush=True)
        result = subprocess.run(
            ["bash", str(source / "regtest/scripts/run")],
            cwd=fixture,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
            timeout=120,
        )
        report_path = fixture / "report.txt"
        report = report_path.read_text() if report_path.exists() else ""
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
                f"Installed regression {name} failed or skipped:\n"
                + report
                + result.stdout
                + result.stderr
            )
        print(f"Installed numerical regression passed: {name}", flush=True)
