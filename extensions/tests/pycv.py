#!/usr/bin/env python3
"""Validate the installed pycv module (toolchain python) against the main PLUMED package."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path, required=True)
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
    PLUMED_KERNEL=str(next((prefix / "lib").glob("libplumedKernel.*"))),
)
subprocess.run([program, "modules", "-q", "--has", "pycv"], env=environment, check=True)
subprocess.run(
    [program, "modules", "-q", "--has-feature", "python_embed"],
    env=environment,
    check=True,
)
assert (
    subprocess.run(
        [program, "modules", "-q", "--has", "pycv"],
        env=dict(environment, PLUMED_AUTOLOAD="0"),
        check=False,
    ).returncode
    == 1
)
subprocess.run(
    [program, "--load-module=pycv", "modules", "-q", "--has", "pycv"],
    env=dict(environment, PLUMED_AUTOLOAD="0"),
    check=True,
)
listed = subprocess.check_output(
    [program, "modules", "--list"], env=environment, text=True
)
assert "pycv_runtime" not in listed.split(), "Private Python helper exposed as a module"
assert not (prefix / "lib/plumed/plugins/pycv_runtime.plumed").exists()
with tempfile.TemporaryDirectory(prefix="plumed-pycv-tests-") as scratch:
    scratch = Path(scratch)
    # Existing Python-hosted tests exercise the companion extension, lifecycle,
    # initialization callbacks, periodicity, atom access and returned values.
    tests = scratch / "pythontests"
    shutil.copytree(source / "regtest/pycv/pythontests", tests)
    subprocess.run(
        [str(prefix / "bin/python"), "-m", "pytest", "-q", str(tests)],
        cwd=scratch,
        env=environment,
        check=True,
        timeout=180,
    )
    regtest = scratch / "regtest"
    regtest.mkdir()
    (regtest / "scripts").symlink_to(
        source / "regtest/scripts", target_is_directory=True
    )
    fixtures = ("rt-getPositions", "rt-func-derivative")
    for mode in ("autoload", "explicit"):
        for name in fixtures:
            fixture = regtest / mode / Path(name).name
            shutil.copytree(source / "regtest/pycv" / name, fixture)
            # Retain calculations, hooks and numerical references. Add only
            # the loading directive to test the explicit module interface.
            if mode == "explicit":
                input_path = fixture / "plumed.dat"
                input_path.write_text("LOAD MODULE=pycv\n" + input_path.read_text())
            result = subprocess.run(
                ["bash", str(source / "regtest/scripts/run")],
                cwd=fixture,
                env=dict(
                    environment,
                    PLUMED_AUTOLOAD="0" if mode == "explicit" else "1",
                    PLUMED_PROGRAM_NAME=program,
                    # the package provides pycv: with PLUMED_AUTOLOAD=0 the check of
                    # plumed_modules=pycv would skip the test (it is loaded by the input)
                    PLUMED_REGTEST_FORCE_REQUIREMENTS="yes",
                ),
                check=False,
                capture_output=True,
                text=True,
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
                        "SKIP_ON_CI",
                    )
                )
            ):
                raise AssertionError(
                    f"{mode} {name}:\n{report}\n{result.stdout}\n{result.stderr}"
                )
            print(
                f"Installed pycv numerical regression passed: {mode} {name}", flush=True
            )
