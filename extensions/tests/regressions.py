#!/usr/bin/env python3
"""Run representative existing numerical fixtures against installed packages."""

import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path)
parser.add_argument("--list-packages", action="store_true")
args = parser.parse_args()
source = Path(__file__).resolve().parents[2]
selection = json.loads((source / "conda/plumed-extensions/regtests.json").read_text())
fixtures = selection["tests"]
if not fixtures:
    parser.error("the conda regression selection must not be empty")
packages = sorted({package for case in fixtures for package in case["packages"]})
if args.list_packages:
    print(" ".join(packages))
    raise SystemExit(0)
if args.prefix is None:
    parser.error("--prefix is required to run tests")
prefix = args.prefix.resolve()
installed = {
    json.loads(record.read_text())["name"]
    for record in (prefix / "conda-meta").glob("*.json")
}
missing = set(packages) - installed
if missing:
    parser.error("missing test packages: " + ", ".join(sorted(missing)))
environment = dict(
    os.environ,
    PLUMED_PROGRAM_NAME=str(prefix / "bin/plumed"),
    PLUMED_AUTOLOAD="1",
    PLUMED_NUM_THREADS="2",
    OMP_NUM_THREADS="2",
    PLUMED_REGTEST_FORCE_REQUIREMENTS="yes"
    if selection.get("force_requirements")
    else "no",
)
with tempfile.TemporaryDirectory(prefix="plumed-installed-regressions-") as scratch:
    regtest = Path(scratch) / "regtest"
    regtest.mkdir()
    for name in ("scripts", "trajectories"):
        (regtest / name).symlink_to(source / "regtest" / name, target_is_directory=True)
    for case in fixtures:
        name = case["regtest"]
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("invalid regtest path: " + name)
        fixture = regtest / name
        shutil.copytree(
            source / "regtest" / name,
            fixture,
            ignore=shutil.ignore_patterns("tmp", "report.txt"),
        )
        # Keep original inputs, hooks and references; eligibility forcing is explicit.
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
