#!/usr/bin/env python3
"""Acceptance checks against an installed base and extension prefix."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--base", required=True, type=Path)
parser.add_argument("--extensions", required=True, type=Path)
parser.add_argument("--visibility", type=Path)
args = parser.parse_args()
base = args.base.resolve()
extensions = args.extensions.resolve()
plumed = str(base / "bin/plumed")
registrations = extensions / "lib/plumed/plugins"
environment = dict(
    os.environ, PLUMED_PLUGIN_PATH=str(registrations), PLUMED_AUTOLOAD="0"
)


def run(*arguments, success=True, env=None, error=None):
    result = subprocess.run(
        [plumed, *arguments],
        env=env or environment,
        text=True,
        capture_output=True,
        check=False,
    )
    if (result.returncode == 0) != success:
        raise AssertionError(
            f"{arguments}: exit {result.returncode}\n{result.stdout}\n{result.stderr}"
        )
    if error and error not in result.stderr:
        raise AssertionError(f"Expected {error!r} in {result.stderr!r}")
    return result.stdout


assert "ves" not in run("--list-modules").split()
run("config", "-q", "module", "ves", success=False)
run("--load-module=ves", "config", "-q", "module", "ves")
run("--load-module=pamm", "config", "-q", "module", "pamm", "adjmat")
run("--load-module=pamm", "config", "-q", "module", "ves", success=False)
run("--load-module=bias", "--load-module=bias", "config", "-q", "module", "bias")
run("--load-module=nonexistent", "config", "-q", "module", "nonexistent", success=False)
run("--load-module=../ves", "help", success=False)
automatic = dict(environment, PLUMED_AUTOLOAD="1")
run("config", "-q", "module", "ves", "pamm", "adjmat", env=automatic)
# Named loading has to reach driver-created PLUMED objects, not just CLI lookup.
with tempfile.TemporaryDirectory() as temp:
    directory = Path(temp)
    input_file = directory / "input.dat"
    input_file.write_text("bf: BF_CHEBYSHEV ORDER=3 MINIMUM=-1 MAXIMUM=1\n")
    run("--load-module=ves", "driver", "--noatoms", "--plumed", str(input_file))
    # A second registration directory cannot silently replace the same module.
    conflict = directory / "conflicts"
    conflict.mkdir()
    manifest = (registrations / "ves.plumed").read_text()
    library_name = next(
        line.split(" ", 1)[1]
        for line in manifest.splitlines()
        if line.startswith("library ")
    )
    real_library = (registrations / library_name).resolve()
    duplicate = directory / real_library.name
    shutil.copy2(real_library, duplicate)
    altered = (
        "\n".join(
            "library " + str(duplicate) if line.startswith("library ") else line
            for line in manifest.splitlines()
        )
        + "\n"
    )
    (conflict / "ves.plumed").write_text(altered)
    conflicting = dict(automatic, PLUMED_PLUGIN_PATH=f"{registrations}:{conflict}")
    run("config", "module", "ves", success=False, env=conflicting)
    selected = dict(conflicting, PLUMED_MODULE_OVERRIDE=f"ves={real_library}")
    run("config", "-q", "module", "ves", env=selected)
    # Build mismatch is rejected before loading the library.
    mismatch = directory / "mismatch"
    mismatch.mkdir()
    (mismatch / "ves.plumed").write_text(
        "\n".join(
            "build incompatible"
            if line.startswith("build ")
            else "module mismatch_ves"
            if line.startswith("module ")
            else line
            for line in altered.splitlines()
        )
        + "\n"
    )
    mismatched = dict(automatic, PLUMED_PLUGIN_PATH=str(mismatch))
    run(
        "config",
        "module",
        "ves",
        success=False,
        env=mismatched,
        error="different PLUMED base",
    )
    # Autoload opt-out must recover even with a bad installed registration.
    run("config", "-q", "module", "bias", env=dict(mismatched, PLUMED_AUTOLOAD="0"))
    absent = dict(automatic, PLUMED_PLUGIN_PATH=str(directory / "absent"))
    run("config", "-q", "module", "bias", env=absent)
if args.visibility:
    subprocess.run(
        [str(args.visibility.resolve())],
        env=dict(
            environment,
            PLUMED_KERNEL=str(next((base / "lib").glob("libplumedKernel.*"))),
        ),
        check=True,
    )
print("Installed extension acceptance checks passed")
